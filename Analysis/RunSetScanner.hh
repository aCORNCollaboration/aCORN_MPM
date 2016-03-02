/// \file RunSetScanner.hh Base class for scanning through files specified by run number
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef RUNSETSCANNER_HH
#define RUNSETSCANNER_HH

#include "TChainScanner.hh"
#include "TCumulativeMap.hh"
#include "Enums.hh"

#include <TH1.h>

#include <map>
#include <cassert>

/// Class for scanning through files specified by run number
class RunSetScanner: public TChainScanner {
public:
    /// constructor
    RunSetScanner(const std::string& treeName): TChainScanner(treeName) { }
    
    /// add run to data (return whether successful)
    virtual bool addRun(RunID rn);
    /// add list of runs to data; return number successfully added
    unsigned int addRuns(const vector<RunID>& rns);
    /// return path to run .root file
    virtual string locateRun(RunID) const { assert(false); return ""; }
    
    /// speedload, keeping track of currently loaded run number
    void speedload(unsigned int e) override;
    /// update calibrations when next TTree loaded
    void nextTreeLoaded() override;
    ///  subclass this for calibrations after loading event
    virtual void calibrate() { }
    /// subclass this for routines when new run is loaded
    virtual void loadNewRun(RunID) { }
    /// get run ID of current event
    virtual RunID getRun() const { return evtRun; }
    /// check whether this is simulated data
    virtual bool isSimulated() const { return false; }
    
    /// print info about this scanner
    virtual void display();
    /// print info about current event
    virtual void displayEvt() const { }
    
    RunID evtRun;                               ///< run ID for current event
    double totalTime = 0;                       ///< combined length of runs in seconds
    TCumulativeMap<RunNum,Double_t>  runTimes;  ///< times for each run loaded
    TCumulativeMap<RunNum,Double_t>  runCounts; ///< total event counts for each run loaded
    
protected:
    
    vector<RunID> runlist;              ///< list of loaded runs
    /// at run load time, figure out run total time
    virtual double _getRunTime(RunID) { return 0; }
};


#endif
