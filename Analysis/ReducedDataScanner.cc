#include "ReducedDataScanner.hh"

void ReducedDataScanner::setReadpoints(TTree* T) {
	T->Branch("nP", &nP, "nP/I");
	T->Branch("nE", &nE, "nE/I");
	T->Branch("V", &V, "V/I");
	
	T->Branch("T_p", &T_p, "T_p/L");
	T->Branch("T_e2p", &T_e2p, "T_e2p/I");
	
	T->Branch("E_p", &E_p, "E_p/I");
	T->Branch("E_e", &E_e, "E_e/I");
	
	T->Branch("DetFired", &DetFired, "DetFired/i");
	T->Branch("DetPiled", &DetPiled, "DetPiled/i");
	T->Branch("flags", &flags, "flags/I");
	T->Branch("Idx", &Idx, "Idx/I");
	
	T->Branch("E_PMT", E_PMT, "E_PMT[32]/S");
	T->Branch("T_PMT_median", &T_PMT_median, "T_PMT_median/B");
	T->Branch("Max_PMT", &Max_PMT, "MaxPMT/B");
	T->Branch("E_Max_PMT", &E_Max_PMT, "E_Max_PMT/S");
}
