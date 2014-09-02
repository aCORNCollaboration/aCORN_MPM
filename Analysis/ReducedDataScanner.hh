#ifndef REDUCEDDATASCANNER_HH
#define REDUCEDDATASCANNER_HH

#include "BaseDataScanner.hh"

/// class for assembling and scanning TChain of "Reduced" data
class ReducedDataScanner: public BaseDataScanner {
public:
    /// constructor
    ReducedDataScanner(bool fp): BaseDataScanner("RedEvt",fp) {}
    
    /// return path to run .root file
    virtual std::string locateRun(RunID r);
    
    /// set TChain branch data readpoints
    virtual void setReadpoints(TTree* T);
    
    /// set TTree write points
    void setWritepoints(TTree* T);
};

#endif
