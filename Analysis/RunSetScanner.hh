#ifndef RUNSETSCANNER_HH
#define RUNSETSCANNER_HH

#include "TChainScanner.hh"
#include "Enums.hh"
#include "AcornCalibrator.hh"

#include "QFile.hh"
#include "TagCounter.hh"

#include <TH1.h>

#include <map>
#include <cassert>

/// Class for scanning through files specified by run number
class RunSetScanner: public TChainScanner {
public:
    /// constructor
    RunSetScanner(const std::string& treeName);
    
    /// destructor
    virtual ~RunSetScanner();
    
    /// add run to data (return whether successful)
    virtual bool addRun(RunID rn);
    /// add list of runs to data; return number successfully added
    unsigned int addRuns(const vector<RunID>& rns);
    /// return path to run .root file
    virtual string locateRun(RunID) { assert(false); return ""; }
    
    /// speedload, keeping track of currently loaded run number
    virtual void speedload(unsigned int e, bool loadBaskets = true);
    ///  subclass this for calibrations after loading event
    virtual void calibrate() { }
    /// subclass this for routines when new run is loaded
    virtual void loadNewRun(RunID) {}
    /// get run ID of current event
    virtual RunID getRun() const { return evtRun; }
    /// check whether this is simulated data
    virtual bool isSimulated() const { return false; }
    
    /// print info about this scanner
    virtual void display();
    
    RunID evtRun;                       ///< run ID for current event
    double totalTime;                   ///< combined length of runs in seconds
    TagCounter<RunID>  runTimes;        ///< times for each run loaded
    TagCounter<RunID>  runCounts;       ///< total event counts for each run loaded
    
protected:
    
    vector<RunID> runlist;                 ///< list of loaded runs
    /// at run load time, figure out run total time
    virtual double _getRunTime(RunID) { return 0; }
};


#endif
