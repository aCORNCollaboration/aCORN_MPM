/// \file ReducedDataScanner.hh Assembles and scans a TChain of "Reduced" data
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef REDUCEDDATASCANNER_HH
#define REDUCEDDATASCANNER_HH

#include "BaseDataScanner.hh"

/// Assembles and scans a TChain of "Reduced" data
class ReducedDataScanner: public BaseDataScanner {
public:
    /// constructor
    ReducedDataScanner(bool fp): BaseDataScanner("RedEvt",fp) {}
    
    /// return path to run .root file
    string locateRun(RunID r) const override;
    
    /// set TChain branch data readpoints
    void setReadpoints(TTree* T) override;
    
    /// set TTree write points
    void setWritepoints(TTree* T);
    
    /// locate runs for specified series on disk
    vector<RunID> findSeriesRuns(int s) const;
};

#endif
