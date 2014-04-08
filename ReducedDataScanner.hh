#ifndef REDUCEDDATASCANNER_HH
#define REDUCEDDATASCANNER_HH

#include "BaseDataScanner.hh"

/// class for assembling and scanning TChain of "Reduced" data
class ReducedDataScanner: public BaseDataScanner {
public:
	/// constructor
	ReducedDataScanner(): BaseDataScanner("RedEvt") {}
	
	/// find path to processed run .root file
	virtual std::string locateRun(RunNum) { return ""; }
		
	/// set TChain branch data readpoints
	virtual void setReadpoints(TTree* T);
};

#endif
