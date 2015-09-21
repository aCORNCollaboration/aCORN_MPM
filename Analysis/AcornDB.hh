/// \file AcornDB.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef ACORNDB_HH
#define ACORNDB_HH

#include "Enums.hh"
#include "CalPeak.hh"
#include "SQLite_Helper.hh"
#include <TGraphErrors.h>
#include <stdlib.h>
#include <vector>

using std::vector;

/// Generic analysis result number
struct AnaResult {
    RunID start;        ///< start run included
    RunID end;          ///< end run included
    time_t time = 0;    ///< timestamp for result record
    double value = 0;   ///< result value
    double err = 0;     ///< result uncertainty
};

class AcornDB: public SQLite_Helper {
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
    /// get list of runs in database for given wishbone series
    vector<RunID> seriesRuns(RunNum S, DataTier T = GOOD);
    /// get list of series in database for given analysis group
    vector<RunNum> groupSeries(const string& sname);
    
    /// upload PMT gain calibration data
    void loadPMTcal(RunID start, RunID end, int n, double sigPerPE, double sigPerMeV);
    /// get PMT gain calibration parameters
    void getPMTcal(RunID rn, vector<double>& sigPerPE, vector<double>& sigPerMeV);
    /// get energy recalibration curve
    TGraphErrors* getRecal(RunID rn);
    
    /// upload peak data, optionally deleting prior matches
    void uploadPeak(const CalPeak& pk, bool replace = true);
    /// create new named object ID
    sqlite3_int64 createNamed(const string& tp, const string& name, const string& descrip);
    /// upload graph (will try dynamic_cast to TGraphErrors); return graph ID
    sqlite3_int64 uploadGraph(const TGraph* g, const string& name, const string& descrip);
    /// get graph by ID
    TGraphErrors* getGraph(sqlite3_int64 gID);
    /// upload energy recalibration curve
    sqlite3_int64 uploadRecal(const TGraph* g, RunID r0, RunID r1);
    
    /// locate / generate analysis result type ID
    sqlite3_int64 getAnaResultType(const string& name, const string& descrip = "");
    /// upload analysis result by type number
    void uploadAnaResult(sqlite3_int64 type_id, AnaResult R);
    /// upload analysis result by name
    void uploadAnaResult(const string& name, const string& descrip, AnaResult R);
    
protected:
    /// Constructor
    AcornDB();

    static AcornDB* myDB;       ///< static singleton instance
};


#endif
