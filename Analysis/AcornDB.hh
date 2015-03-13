#ifndef ACORNDB_HH
#define ACORNDB_HH

#include "Enums.hh"
#include "CalPeak.hh"
#include "SQLite_Helper.hh"
#include <stdlib.h>
#include <vector>

using std::vector;

class AcornDB: protected SQLite_Helper {
public:
    /// get singleton instance
    static AcornDB& ADB();
    
    enum DataTier {
        UNKNOWN   = 0,
        GOOD      = 1,
        USABLE    = 2,
        UNUSUABLE = 3
    };
    
    /// get PMT sum energy calibration parameters
    void getPMTSumCal(RunID rn, double& slope, double& intercept);
    /// get list of runs for given wishbone series
    vector<RunID> seriesRuns(RunNum S, DataTier T = GOOD);
    
    /// upload PMT gain calibration data
    void loadPMTcal(RunID start, RunID end, int n, double sigPerPE, double sigPerMeV);
    /// get PMT gain calibration parameters
    void getPMTcal(RunID rn, vector<double>& sigPerPE, vector<double>& sigPerMeV);
    
    /// upload peak data, optionally deleting prior matches
    void uploadPeak(const CalPeak& pk, bool replace = true);
    
protected:
    /// Constructor
    AcornDB();
    
    static AcornDB* myDB;       ///< static singleton instance
};


#endif
