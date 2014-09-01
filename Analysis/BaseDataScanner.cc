#include "BaseDataScanner.hh"

BaseDataScanner::BaseDataScanner(const std::string& treeName):
RunSetScanner(treeName), physicsWeight(1.) { }

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
    if(E_p != -1) E_p_0 = E_p;
    T_p *= NS_PER_CLOCK;
    T_e2p *= NS_PER_CLOCK;
    E_recon = currentCal->calEnergy(E_e);
}
