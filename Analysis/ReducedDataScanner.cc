/// \file ReducedDataScanner.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "ReducedDataScanner.hh"
#include "PathUtils.hh"
#include "StringManip.hh"
#include <stdio.h>

void ReducedDataScanner::setReadpoints(TTree* T) {
    SetBranchAddress(T, "nP", &nP);
    if(is4p) SetBranchAddress(T, "nPSig", &nPSig);
    SetBranchAddress(T, "nE", &nE);
    SetBranchAddress(T, "nV", &nV);
    
    SetBranchAddress(T, "T_p", &T_p);
    if(is4p) SetBranchAddress(T, "T_d", &T_d);
    SetBranchAddress(T, "T_e2p", &T_e2p);
    if(is4p) SetBranchAddress(T, "T_e2d", &T_e2d);
    
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
    if(is4p) T->Branch("nPSig", &nPSig, "nPSig/I");
    T->Branch("nE", &nE, "nE/I");
    T->Branch("nV", &nV, "nV/I");
    
    T->Branch("T_p", &T_p, "T_p/L");
    if(is4p) T->Branch("T_d", &T_d, "T_d/L");
    T->Branch("T_e2p", &T_e2p, "T_e2p/I");
    if(is4p) T->Branch("T_e2d", &T_e2d, "T_e2d/I");
        
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

string ReducedDataScanner::locateRun(RunID r) const {
    char fnm[1024];
    sprintf(fnm,"%s/s%ir%04i_red.root", getEnvSafe("ACORN_REDUCED_ROOT").c_str(), r.first, r.second);
    if(fileExists(fnm)) return fnm;
    sprintf(fnm,"%s/s%04ir%04i_red.root", getEnvSafe("ACORN_REDUCED_ROOT").c_str(), r.first, r.second);
    if(fileExists(fnm)) return fnm;
    sprintf(fnm,"%s/s%04ir%04i_rd2.root", getEnvSafe("ACORN_REDUCED_ROOT").c_str(), r.first, r.second);
    if(fileExists(fnm)) return fnm;
    return "";
}

vector<RunID> ReducedDataScanner::findSeriesRuns(int s) const {
    vector<RunID> v;
    unsigned int nmissing = 0;
    RunID rn;
    rn.first = s;
    rn.second = 0;
    while(nmissing < 10) {
        if(locateRun(rn).size()) {
            v.push_back(rn);
            nmissing = 0;
        } else nmissing++;
        rn.second++;
    }
    return v;
}

