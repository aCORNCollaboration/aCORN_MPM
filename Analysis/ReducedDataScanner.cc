#include "ReducedDataScanner.hh"
#include "PathUtils.hh"
#include "strutils.hh"
#include <stdio.h>

void ReducedDataScanner::setReadpoints(TTree* T) {
    SetBranchAddress(T, "nP", &nP);
    SetBranchAddress(T, "nE", &nE);
    SetBranchAddress(T, "V", &V);
    
    SetBranchAddress(T, "T_p", &T_p);
    SetBranchAddress(T, "T_e2p", &T_e2p);
    
    SetBranchAddress(T, "E_p", &E_p);
    SetBranchAddress(T, "E_e", &E_e);
    
    SetBranchAddress(T, "DetFired", &DetFired);
    SetBranchAddress(T, "DetPiled", &DetPiled);
    SetBranchAddress(T, "flags", &flags);
    SetBranchAddress(T, "Idx", &Idx);
    
    SetBranchAddress(T, "E_PMT", E_PMT);
    SetBranchAddress(T, "T_PMT_median", &T_PMT_median);
    SetBranchAddress(T, "Max_PMT", &Max_PMT);
    SetBranchAddress(T, "E_Max_PMT", &E_Max_PMT);
}

void ReducedDataScanner::setWritepoints(TTree* T) {
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

std::string ReducedDataScanner::locateRun(RunID r) {
    char fnm[1024];
    sprintf(fnm,"%s/s%04ir%04i_red.root", getEnvSafe("ACORN_REDUCED_ROOT").c_str(), r.first, r.second);
    //printf("Looking for '%s'...\n",fnm);
    if(fileExists(fnm)) return fnm;
    sprintf(fnm,"%s/s%04ir%04i_rd2.root", getEnvSafe("ACORN_REDUCED_ROOT").c_str(), r.first, r.second);
    if(fileExists(fnm)) return fnm;
    return "";
}