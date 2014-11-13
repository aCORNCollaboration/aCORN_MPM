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
    
    /// Individual-PMT optimal sum calibration
    double calEnergy(const Short_t* ADC) const;
    
protected:
    RunID rn;           ///< Run ID being calibrated
    
    double slope;       ///< simple calibration slope
    double intercept;   ///< simple calibration intercept
    
    vector<double> sigPerPE;       ///< energy resolution data
    vector<double> sigPerMeV;      ///< gain calibration data
    double PEPerMeV;                    ///< calibrated total PE per MeV
};

#endif
