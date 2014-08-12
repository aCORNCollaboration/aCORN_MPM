#ifndef RUNMETADATA_HH
#define RUNMETADATA_HH

#include "Enums.hh"
#include <map>
#include <vector>

using std::map;
using std::vector;

/// Metadata for aCORN run
class RunMetadata {
public:
    /// Constructor
    RunMetadata(): tier(UNKNOWN), field(FIELD_UNKNOWN), slope(0), intercept(0) {}
    
    RunID run;  ///< run identifier
    
    enum DataTier {
        UNKNOWN   = 0,
        GOOD      = 1,
        USABLE    = 2,
        UNUSUABLE = 3
    } tier;     ///< quality-of-data tier

    enum Field {
        FIELD_DOWN    = -1,
        FIELD_UNKNOWN = 0,
        FIELD_UP      = 1
    } field;    ///< field direction
    
    /// Calibrate PMT sum to energy
    double calibPMT(double S) const { return slope*S + intercept; }
    
    double slope;       ///< calibration slope
    double intercept;   ///< calibration intercept
};

#include "PathUtils.hh"
#include <string>

/// Metadata file loader
class MetadataDB {
public:
    
    /// get metadata for run
    RunMetadata getRun(RunID r) const;
    
    /// get runlist by tier for series
    vector<RunID> seriesRuns(RunNum S, RunMetadata::DataTier T = RunMetadata::GOOD) const;
    
    static MetadataDB MDB;      /// singleton instance
    
protected:
    
    /// Constructor, optionally loading CSV file
    MetadataDB(const std::string& fname = getEnvSafe("ACORN_METADATA"));
    
    map<RunID, RunMetadata> md;                 ///< available metadata
    map<RunNum, vector<RunNum>> series;         ///< listing by series
};

#endif
