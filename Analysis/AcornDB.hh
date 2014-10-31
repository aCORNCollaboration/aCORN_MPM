#ifndef ACORNDB_HH
#define ACORNDB_HH

#include "Enums.hh"
#include "sqlite3.h"
#include <stdlib.h>
#include <vector>

using std::vector;

class AcornDB {
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
    
    /// upload PMT calibration data
    void loadPMTcal(RunID start, RunID end, int n, double sigPerPE, double sigPerMeV);
    
protected:
    /// Constructor
    AcornDB();
    /// Destructor
    ~AcornDB();
    
    /// set up query for use
    int setQuery(const char* qry, sqlite3_stmt*& stmt);
    
    sqlite3* db = NULL;         ///< database connection
    static AcornDB* myDB;       ///< static singleton instance
    vector<sqlite3_stmt*> statements;    ///< prepared statements awaiting deletion
};


#endif
