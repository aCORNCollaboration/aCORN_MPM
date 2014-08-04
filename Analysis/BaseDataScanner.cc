#include "BaseDataScanner.hh"

BaseDataScanner::BaseDataScanner(const std::string& treeName):
RunSetScanner(treeName) { }

void BaseDataScanner::makeFlags() {
    flags = TriggerCategory(0);
    if(500 < E_p && E_p < 2000) flags = TriggerCategory(flags | TCAT_PROTON);
    if(7000 < E_p && E_p < 9000) flags = TriggerCategory(flags | TCAT_PULSER);
    if(250 < T_e2p && T_e2p < 450) flags = TriggerCategory(flags | TCAT_TOF_WB);
}
