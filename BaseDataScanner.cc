#include "BaseDataScanner.hh"
#include "PathUtils.hh"
#include "SMExcept.hh"
#include <stdio.h>
#include <stdlib.h>


BaseDataScanner::BaseDataScanner(const std::string& treeName):
RunSetScanner(treeName) { }

Stringmap BaseDataScanner::evtInfo() {
	Stringmap m;
	return m;
}

void BaseDataScanner::makeFlags() {
	flags = 0;
	if(500 < E_p && E_p < 2000) flags |= TCAT_PROTON;
	if(7000 < E_p && E_p < 9000) flags |= TCAT_PULSER;
	if(250 < T_e2p && T_e2p < 450) flags |= TCAT_TOF_WB;
}
