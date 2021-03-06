/// \file BaseDataScanner.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "BaseDataScanner.hh"
#include <iostream>
#include <stdio.h>

BaseDataScanner::BaseDataScanner(const std::string& treeName, bool fp):
RunSetScanner(treeName), is4p(fp), physicsWeight(1.) { 
    
    // 0b0101000000000111111111100001111,
    // 0b1010111111111000000000011110000
    //const UInt_t mod_mask[N_MODULES] = { 671350543, 1476133104 };
    // 0b0101000000000111111111100000000,
    // 0b1010111111111000000000000000000
    const UInt_t mod_mask[N_MODULES] = { 671350528, 1476132864 };
    for(size_t i=0; i<NCH_MAX; i++) module_map[i] = N_MODULES;
    for(unsigned int m=0; m<N_MODULES; m++) {
        for(unsigned int c=0; c<NCH_MAX; c++)
            if(mod_mask[m] & 1<<c) module_map[c] = m;
    }
    
}

BaseDataScanner::~BaseDataScanner() {
    for(auto kv: cals) delete kv.second;
}

void BaseDataScanner::makeFlags() {
    flags = TriggerCategory(0);
    if(500 < E_p && E_p < 2000) flags = TriggerCategory(flags | TCAT_PROTON);
    if(7000 < E_p && E_p < 9000) flags = TriggerCategory(flags | TCAT_PULSER);
    if(250 < T_e2p && T_e2p < 450) flags = TriggerCategory(flags | TCAT_TOF_WB);
}

void BaseDataScanner::loadNewRun(RunID rn) {
    if(evtRun < rn) T_0 += runTimes[toRunNum(evtRun)] * 1e9;
    auto it = cals.find(rn);
    if(it == cals.end()) {
        currentCal = cals[rn] = new AcornCalibrator(rn);
    } else {
        currentCal = it->second;
    }
}

double BaseDataScanner::_getRunTime(RunID) {
    int nback = 1;
    T_p = 0;
    while(T_p <= 0) speedload(nEvents-(nback++));
    return T_p > 0 ? T_p/1.e9 : 144.; // TODO correct this!
}

void BaseDataScanner::calibrate() {
    if(E_p > 0) E_p_0 = E_p;
    T_p *= NS_PER_CLOCK;
    T_d *= NS_PER_CLOCK;
    T_e2p *= NS_PER_CLOCK;
    T_e2d *= NS_PER_CLOCK;
    
    assert(currentCal);
    //E_recon = currentCal->calPMTSum(E_e); // TODO auto switch wishbone/sourcecal
    E_recon = currentCal->calWishboneSum(E_e);
    
    for(unsigned int i=0; i<N_E_PMT; i++) L_PMT[i] = E_PMT[N_V_PMT+i];
    Pos.calcPos(L_PMT);
    
    nFiredModule();
    modDropoutEvt = !(nFiredMod[0]*nFiredMod[1]) && (nFiredMod[0]+nFiredMod[1]>=6);
}

void BaseDataScanner::nFiredModule() {
    for(unsigned int m=0; m<=N_MODULES; m++) nFiredMod[m] = 0;
    for(unsigned int i=0; i<NCH_MAX; i++) if(DetFired & (1<<i)) nFiredMod[module_map[i]]++;
}

void BaseDataScanner::displayEvt() const {
    std::cout << "E_p = " << E_p << "\tE_p_0 = " << E_p_0 << "\tE_e = " << E_e << "\tT_p = " << T_p << "\n";
}
