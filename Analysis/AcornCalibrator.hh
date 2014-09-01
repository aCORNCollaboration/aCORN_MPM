#ifndef ACORNCALIBRATOR_HH
#define ACORNCALIBRATOR_HH

#include "Enums.hh"
#include "RunMetadata.hh"
#include "ManualInfo.hh"

/// Calibrator and associated information class for aCORN data
class AcornCalibrator {
public:
    /// Constructor
    AcornCalibrator(RunID r);
    
    /// Calibrate electron energy
    double calEnergy(double ADC) const;
    
protected:
    RunID rn;           ///< Run ID being calibrated
    RunMetadata md;     ///< metadata file entry for run
};

#endif
