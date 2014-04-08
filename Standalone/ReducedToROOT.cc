#include <RTypes.h>
#include <TTree.h>
#include <TFile.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <cassert>

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
//	5-digit time offset between electron and proton energy signal (eTIme = P_Time - EPTime)
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
	const unsigned int NCH_MAX = 32;
	Long_t P_Time;			// proton time, 10ns
	Int_t PEngy;			// proton energy
	Int_t nP = 1;			// number of protons (always 1 in this format)
	Int_t PmETime;			// P_Time - EPTime
	Int_t E_Enrgy;			// Electron energy
	UInt_t nE;				// number of electron PMTs
	Int_t V;				// number of veto PMTs
	UInt_t DetFired;		// detector firing bit map 0-30
	UInt_t DetPiled;		// pileup flags [e- dead zone][det. 0-30]
	Int_t Idx;				// proton index
	UInt_t dtl[NCH_MAX];	// Electron PMT "detail" code
	TTree* T  = new TTree("RedEvt","aCORN reduced event data");
	T->Branch("P_Time", &P_Time, "P_Time/L");
	T->Branch("PEngy", &PEngy, "PEngy/I");
	T->Branch("nP", &nP, "nP/I");
	T->Branch("PmETime", &PmETime, "PmETime/I");
	T->Branch("E_Enrgy", &E_Enrgy, "E_Enrgy/I");
	T->Branch("nE", &nE, "nE/I");
	T->Branch("V", &V, "V/I");
	T->Branch("DetFired", &DetFired, "DetFired/i");
	T->Branch("DetPiles", &DetPiled, "DetPiled/i");
	T->Branch("Idx", &Idx, "Idx/I");
	T->Branch("dtl", dtl, "dtl[32]/i");
	
	while(!inf.fail()) {
		inf >> std::dec >> std::skipws >>  P_Time >> PEngy >> PmETime >> E_Enrgy >> nE >> V
			>> std::hex >> DetFired >> DetPiled >> Idx;
		
		if(inf.fail()) break;
		
		//std::cout << P_Time << "\t" << PEngy << "\t" << nE;
		assert(nE+V < NCH_MAX);
		for(unsigned int i=0; i<nE+V; i++) {
			inf >> std::hex >> dtl[i];
			//std::cout << "\t" << dtl[i];
		}
		//std::cout << "\n";
		
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
#####  P_Time   PEngy E-PTm E_Enrgy nE V DetFired DetPiled --Index-  EPMT codes
000000000066474   320     0       0  0 0 28000000 00000000 00000000

000000007529767 14137   600    5177 11 0 2FF1C800 00000000 00000004 700C90 780C42 D00B12 B03712 581613 B80C43 C81C04 C01854 A828A4 A03E14 802B55
*/

/*
TChain T("RedEvt")
T.Add("*.root") 
# proton spectrum; main peak 500--2000; pulser at 7000--9000
T->Draw("PEngy","PEngy<10000")

# electron energy? so badly uncalibrated?
T->Draw("E_Enrgy","PEngy>500 && PEngy<2000 && V==0 && E_Enrgy > 10 && E_Enrgy < 2000")
T->Draw("E_Enrgy","PEngy>500 && PEngy<2000 && V>0 && E_Enrgy > 10 && E_Enrgy < 2000","Same")

# proton event rate
T->Draw("P_Time*(1e-9)","PEngy>500 && PEngy<2000 && V==0")

# Wishbone!!
T->Draw("PmETime*0.01:E_Enrgy","PEngy>500 && PEngy<2000 && V==0 && E_Enrgy > 4000 && E_Enrgy < 35000 && PmETime*0.01 > 2.5 && PmETime*0.01 < 4.5","Col")

TH2F foo("foo","aCORN wishbone",50,4000,35000,50,2.5,4.5)
T->Draw("PmETime*0.01:E_Enrgy >> foo","PEngy>500 && PEngy<2000 && V==0")

gStyle->SetOptStat("")
foo.GetXaxis()->SetTitle("e^{-} energy [uncalibrated ADC]")
foo.GetYaxis()->SetTitle("proton TOF [#mus]")
foo.Draw("Col")
*/
