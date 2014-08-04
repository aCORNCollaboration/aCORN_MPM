#include <RTypes.h>
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

// Reduce14 fields:
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
//	"Detail" codes
//


int main(int argc, char** argv) {
    if(argc<3) {
        printf("Use: ./ReducedToROOT [reduced.txt filename] [output.root filename]\n");
        return 1;
    }
    
    // input file
    std::ifstream inf;
    inf.open(argv[1]);
    
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
    ReducedDataScanner R;
    R.setWritepoints(T);
    
    while(!inf.fail()) {
        inf >> std::dec >> std::skipws >>  R.T_p >> R.E_p >> R.T_e2p >> R.E_e >> R.nE >> R.V
        >> std::hex >> R.DetFired >> R.DetPiled >> R.Idx;
        
        if(inf.fail()) break;
        
        for(unsigned int i=0; i<NCH_MAX; i++) R.E_PMT[i] = R.T_PMT[i] = 0;
        
        //std::cout << T_p << "\t" << E_p << "\t" << nE;
        assert(R.nE+R.V < NCH_MAX);
        UInt_t dtl[NCH_MAX];	// Electron PMT "detail" code
        
        std::vector<Char_t> tps;
        R.Max_PMT = -1;
        R.E_Max_PMT = 0;
        for(unsigned int i=0; i<R.nE+R.V; i++) {
            inf >> std::hex >> dtl[i];
            //std::cout << "\t" << dtl[i];
            unsigned int chn = (dtl[i]>>19);
            assert(chn < NCH_MAX);
            R.E_PMT[chn] = (dtl[i]>>4) & ((1<<15)-1);
            R.T_PMT[chn] = dtl[i] & ((1<<4)-1);
            tps.push_back(R.T_PMT[chn]);
            if(R.E_Max_PMT < R.E_PMT[chn]) { R.Max_PMT = chn; R.E_Max_PMT = R.E_PMT[chn]; }
        }
        //std::cout << "\n";
        std::sort(tps.begin(),tps.end());
        R.T_PMT_median = tps.size()?tps[tps.size()/2]:0;
        R.makeFlags();
        
        // make sure we're at the end of the line (note '\r' endings...)
        inf.getline(lbuf,2,'\r');
        assert(!inf.fail());
        
        T->Fill();
    }
    
    std::cout << "Loaded " << T->GetEntries() << " events.\n";
    
    T->Write();
    outf->Close();
    return 0;
}

/*
#####  T_p   E_p E-PTm E_e nE V DetFired DetPiled --Index-  EPMT codes
000000000066474   320     0       0  0 0 28000000 00000000 00000000

000000007529767 14137   600    5177 11 0 2FF1C800 00000000 00000004 700C90 780C42 D00B12 B03712 581613 B80C43 C81C04 C01854 A828A4 A03E14 802B55

The top 5 bits are the detector number.
The next 11 or15 bits are the individual PMT energy. The original PIXIE value was a 15-bit number but early versions of the reduced format used only the 11 most significant bits to save some space.
The final 4 bits are the time offset (in 10nS units) of this individual PMT relative to the last PMT event in the electron. The actual firing time of this individual PMT can be found from
PMT Time = Proton Time - E-P Delay (from field3) - time offset.

*/

/*
TChain T("RedEvt")
T.Add("*.root") 
# proton spectrum; main peak 500--2000; pulser at 7000--9000
T->Draw("E_p","E_p<10000")

# electron energy? so badly uncalibrated?
T->Draw("E_e","E_p>500 && E_p<2000 && V==0 && E_e > 10 && E_e < 2000")
T->Draw("E_e","E_p>500 && E_p<2000 && V>0 && E_e > 10 && E_e < 2000","Same")

# proton event rate
T->Draw("T_p*(1e-9)","E_p>500 && E_p<2000 && V==0")

# Wishbone!!
T->Draw("T_e2p*0.01:E_e","E_p>500 && E_p<2000 && V==0 && E_e > 4000 && E_e < 35000 && T_e2p*0.01 > 2.5 && T_e2p*0.01 < 4.5","Col")

TH2F foo("foo","aCORN wishbone",50,4000,35000,50,2.5,4.5)
T->Draw("T_e2p*0.01:E_e >> foo","E_p>500 && E_p<2000 && V==0")

gStyle->SetOptStat("")
foo.GetXaxis()->SetTitle("e^{-} energy [uncalibrated ADC]")
foo.GetYaxis()->SetTitle("proton TOF [#mus]")
foo.Draw("Col")

# multiplicity
T->Draw("nE:E_e","E_p>500 && E_p<2000 && V==0 && E_e > 1000 && E_e < 35000 && T_e2p*0.01 > 2.5 && T_e2p*0.01 < 4.5","Col")

# sum is biased!
T.Draw("E_e : (E_PMT[0]+E_PMT[1]+E_PMT[2]+E_PMT[3] + E_PMT[4]+E_PMT[5]+E_PMT[6]+E_PMT[7] + E_PMT[8]+E_PMT[9]+E_PMT[10]+E_PMT[11] + E_PMT[12]+E_PMT[13]+E_PMT[14]+E_PMT[15] + E_PMT[16]+E_PMT[17]+E_PMT[18] + E_PMT[19]+E_PMT[20]+E_PMT[21]+E_PMT[22] + E_PMT[23]+E_PMT[24]+E_PMT[25]+E_PMT[26])*1.0 - E_e","E_e>0 && V==0","Col")
*/
