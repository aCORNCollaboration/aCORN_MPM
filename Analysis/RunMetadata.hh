#ifndef RUNMETADATA_HH
#define RUNMETADATA_HH

#include "Enums.hh"

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
    double calibPMT(double S) { return slope*S + intercept; }
    
    double slope;       ///< calibration slope
    double intercept;   ///< calibration intercept
};

#include "PathUtils.hh"
#include <map>
#include <string>

/// Metadata file loader
class MetadataDB {
public:
    /// Constructor, optionally loading CSV file
    MetadataDB(const std::string& fname = getEnvSafe("ACORN_METADATA"));
    
    /// get metadata for run
    RunMetadata getRun(RunID r) const;
    
protected:
    std::map<RunID, RunMetadata> md;                    ///< available metadata
    std::map<RunNum, std::vector<RunNum>> series;       ///< listing by series
};

#endif
