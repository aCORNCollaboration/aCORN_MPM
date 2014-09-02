#include "BaseDataScanner.hh"

BaseDataScanner::BaseDataScanner(const std::string& treeName, bool fp):
RunSetScanner(treeName), is4p(fp), physicsWeight(1.) { }

void BaseDataScanner::makeFlags() {
    flags = TriggerCategory(0);
    if(500 < E_p && E_p < 2000) flags = TriggerCategory(flags | TCAT_PROTON);
    if(7000 < E_p && E_p < 9000) flags = TriggerCategory(flags | TCAT_PULSER);
    if(250 < T_e2p && T_e2p < 450) flags = TriggerCategory(flags | TCAT_TOF_WB);
}

void BaseDataScanner::loadNewRun(RunID rn) {
    auto it = cals.find(rn);
    if(it == cals.end())
        currentCal = cals[rn] = new AcornCalibrator(rn);
    else
        currentCal = it->second;
}

void BaseDataScanner::calibrate() {
    if(E_p > 0) E_p_0 = E_p;
    T_p *= NS_PER_CLOCK;
    T_d *= NS_PER_CLOCK;
    T_e2p *= NS_PER_CLOCK;
    T_e2d *= NS_PER_CLOCK;
    E_recon = currentCal->calEnergy(E_e);
    nFiredModule();
    modDropoutEvt = !(nFiredMod[0]*nFiredMod[1]) && (nFiredMod[0]+nFiredMod[1]>=6);
}

void BaseDataScanner::nFiredModule() {
    // 0b0101000000000111111111100001111,
    // 0b1010111111111000000000011110000
    //const UInt_t mod_mask[N_MODULES] = { 671350543, 1476133104 };
    // 0b0101000000000111111111100000000,
    // 0b1010111111111000000000000000000
    const UInt_t mod_mask[N_MODULES] = { 671350528, 1476132864 };
    for(unsigned int m=0; m<N_MODULES; m++) {
        nFiredMod[m] = 0;
        for(unsigned int i=0; i<NCH_MAX; i++)
            nFiredMod[m] += ((mod_mask[m] & (1<<i)) && (DetFired & (1<<i)));
    }
}
