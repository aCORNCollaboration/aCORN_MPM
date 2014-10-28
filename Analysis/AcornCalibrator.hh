#ifndef ACORNCALIBRATOR_HH
#define ACORNCALIBRATOR_HH

#include "Enums.hh"
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
    
    double slope;       ///< simple calibration slope
    double intercept;   ///< simple calibration intercept
};

#endif
