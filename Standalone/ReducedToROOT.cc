#include <Rtypes.h>
#include <TTree.h>
#include <TFile.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <cassert>
#include <vector>
#include <algorithm>
#include "ReducedDataScanner.hh"

// ./ReducedToROOT ~/Documents/aCORN_data/s1331pt1/s1331r0000_red.txt ~/Documents/aCORN_data/s1331r0000_red.root


// .red.txt 2p format:
// # #*###  T_p   E_p E-PTm E_e nE V DetFired DetPiled --Index-  EPMT codes
// 000000000066474   320     0       0  0 0 28000000 00000000 00000000
// 000000007529767 14137   600    5177 11 0 2FF1C800 00000000 00000004 700C90 780C42 D00B12 B03712 581613 B80C43 C81C04 C01854 A828A4 A03E14 802B55
// 
// The top 5 bits are the detector number.
// The next 11 or15 bits are the individual PMT energy. The original PIXIE value was a 15-bit number but early versions of the reduced format used only the 11 most significant bits to save some space.
// The final 4 bits are the time offset (in 10nS units) of this individual PMT relative to the last PMT event in the electron. The actual firing time of this individual PMT can be found from
// PMT Time = Proton Time - E-P Delay (from field3) - time offset.
 
 
// Reduce14 (.rd2) fields:
//	15-digit energy arrival time, or -1 if missing
//	15-digit discriminator arrival time, or -1 if missing
//	5-digit proton energy
//		The proton energy could be missing for three reasons which are distinguished by the energy value
//		PEnergy=0 the proton was piled up
//		PEnergy=-1 flag for second and more electrons associated with a single proton.
//		PEnergy=-2 no energy signal detected, event had only discriminator signals.
//	1-digit number of proton signals
//	5-digit time offset between electron and proton energy signal (eTIme = T_p - EPTime)
//	5-digit time offset between electron and proton discriminator signal (eTime = D_Time - EPTime)
//	5-digit electron energy
//	2 digit number of electron PMTs
//	1 digit number of veto PMTs
//
//	8-hex-digit bit map of which detectors (0-30) fired.
//	8-hex-digit bitmap, 31 pile up flags (0-30) and top bit marks events in the dead zone after an electron.
//	8-hex digit proton index.
//
//      The proton details for energy channels are constructed exactly the same as electron details,
//      5 bits of channel #, 15 bits of energy, 4 bits of time offset.
//      The proton details for discriminator channels are constructed differently,
//      5 bits of channel #, 15 bits of time offset, 4 bits of 0
//      The discriminator details are different because there is no energy and the time
//      offset falls outside the 150nS window possible with a 4-bit time offset so I used
//      the energy field to code the time offset.
//      UNLIKE the electron time offsets, the proton time offsets are all time AFTER the first
//      signal that started the event.
//
// #####  P_Time          D_Time   PEngy nP EPTm  EDTm E_Enrgy nE V DetFired DetPiled --Index- PMT codes
// -00000000000001 000000000018330    -2 2     0   514     306  2 0 50001400 00000000 00000000 600550 500DD1 F004C0 E004C0 

int main(int argc, char** argv) {
    if(argc<3) {
        printf("Use: ./ReducedToROOT [reduced.txt filename] [output.root filename]\n");
        return 1;
    }
    
    size_t nsz = strlen(argv[1]);
    if(nsz < 9) {
        std::cout << "Badly formatted input filename; aborting!\n";
        return 1;
    }
    bool is4p = (argv[1][nsz-5]=='2');
    
    // input file
    std::ifstream inf;
    inf.open(argv[1]);
    std::cout << "Processing " << (is4p?"4p '":"2p '") << argv[1] << "'...\n";
    
    // scan to "####" block starting data
    char lbuf[1024];
    while(!inf.fail()) {
        inf.getline(lbuf,1023,'\n');
        if(lbuf[0]=='#') break;
    }
    if(lbuf[0] != '#') {
        printf("Data block not found!\n");
        return -1;
    }
    
    // output file
    TFile* outf = new TFile(argv[2],"RECREATE");
    outf->cd();
    
    // data TTree
    TTree* T  = new TTree("RedEvt","aCORN reduced event data");
    ReducedDataScanner R(is4p);
    R.setWritepoints(T);
    R.nPSig = 0;
    Long64_t T_p_max = 0;
    
    int badLines = 0;
    while(inf.good()) {
        int stage = 0;
        if(is4p) {
            while(1) {
                inf >> std::dec >> std::skipws >> R.T_p; stage++; if(!inf.good()) break;
                inf >> std::dec >> std::skipws >> R.T_d; stage++; if(!inf.good()) break;
                inf >> std::dec >> std::skipws >> R.E_p; stage++; if(!inf.good()) break;
                inf >> std::dec >> std::skipws >> R.nPSig; stage++; if(!inf.good()) break;
                inf >> std::dec >> std::skipws >> R.T_e2p; stage++; if(!inf.good()) break;
                inf >> std::dec >> std::skipws >> R.T_e2d; stage++; if(!inf.good()) break;
                inf >> std::dec >> std::skipws >> R.E_e; stage++; if(!inf.good()) break;
                inf >> std::dec >> std::skipws >> R.nE; stage++; if(!inf.good()) break;
                inf >> std::dec >> std::skipws >> R.nV; stage++; if(!inf.good()) break;
                inf >> std::hex >> std::skipws >>  R.DetFired >> R.DetPiled >> R.Idx;
                break;
            }
        } else {
            inf >> std::dec >> std::skipws >>  R.T_p >> R.E_p >> R.T_e2p >> R.E_e >> R.nE >> R.nV
            >> std::hex >> R.DetFired >> R.DetPiled >> R.Idx;
        }
        
        if(!inf.good()) {
            badLines++;
            //std::cout << "Read fail " << inf.rdstate() << " at stage " << stage << "\n";
            inf.clear();
            inf.getline(lbuf, 512, '\r');
            //std::cout << "Bad line '" << lbuf << "'\n";
            continue;
        }
        
        for(unsigned int i=0; i<NCH_MAX; i++) R.E_PMT[i] = R.T_PMT[i] = 0;
        
        //std::cout << T_p << "\t" << E_p << "\t" << nE;
        //std::cout << "nE = " << R.nE << " nV = " << R.nV << "\n";
        assert(R.nE+R.nV < NCH_MAX);
        
        // Electron PMT "detail" codes
        std::vector<Char_t> tps;
        R.Max_PMT = -1;
        R.E_Max_PMT = 0;
        for(unsigned int i=0; i<R.nE+R.nV+R.nPSig; i++) {
            UInt_t dtl;
            inf >> std::hex >> dtl;
            unsigned int chn = (dtl>>19);               // channel number
            if(chn >= NCH_MAX) continue;                // skip proton readouts, for now...
            R.E_PMT[chn] = (dtl>>4) & ((1<<15)-1);      // PMT ADC
            R.T_PMT[chn] = dtl & ((1<<4)-1);            // PMT TDC
            tps.push_back(R.T_PMT[chn]);
            if(R.E_Max_PMT < R.E_PMT[chn]) { R.Max_PMT = chn; R.E_Max_PMT = R.E_PMT[chn]; }
        }
        //std::cout << "\n";
        std::sort(tps.begin(),tps.end());
        R.T_PMT_median = tps.size()?tps[tps.size()/2]:0;
        R.makeFlags();
        
        // make sure we're at the end of the line (note '\r' endings...)
        inf.getline(lbuf,2,'\r');
        if(inf.fail()) {
            std::cout << "'" << argv[1] << "' hit mal-formed line ending!\n";
            break;
        }
        
        if(T_p_max < R.T_p) T_p_max = R.T_p;
        T->Fill();
    }
    
    std::cout << "'" << argv[1] << "' Loaded " << T->GetEntries() << " events; " << badLines << " bad lines; t= " << T_p_max << "\n";
    
    T->Write();
    outf->Close();
    return 0;
}
