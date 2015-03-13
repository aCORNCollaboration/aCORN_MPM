#ifndef ACORNCALIBRATOR_HH
#define ACORNCALIBRATOR_HH

#include "Enums.hh"
#include "CalPeak.hh"
#include "ManualInfo.hh"

/// Calibrator and associated information class for aCORN data
class AcornCalibrator {
public:
    /// Constructor
    AcornCalibrator(RunID r);
    
    /// Calibrate electron energy from PMT sum
    double calPMTSum(double ADC) const { return slope*ADC + intercept; }
    /// Invert PMT sum calibration
    double invcalPMTSum(double E) const { return (E-intercept)/slope; }
    /// calibration derivative
    double dEds(double) const { return slope; }
    /// Invert PMT sum calibration in peak fit
    void invcalPMTSum(CalPeak& pk) const;
    
    /// optimal weighted sum (effective PE units) from ADCs
    double wsum(const Short_t* ADC) const;
    /// wsum to energy
    double calWSum(double w) const { return  1000 * w / PEPerMeV; }
    /// invert energy to wsum
    double invcalWSum(double E) const { return E/1000 * PEPerMeV; }
    /// calibration derivative
    double dEdw(double) const { return 1000 / PEPerMeV; }
    /// Invert PMT weight sum calibration in peak fit
    void invcalWSum(CalPeak& pk) const;
    
    /// Individual-PMT optimal sum calibration
    double calEnergy(const Short_t* ADC) const { return calWSum(wsum(ADC)); }
    
    /// whether PMT is in outer ring
    static bool isOuter(size_t i);
    
protected:
    RunID rn;                   ///< Run ID being calibrated
    
    double slope;               ///< simple calibration slope
    double intercept;           ///< simple calibration intercept
    
    vector<double> sigPerPE;    ///< energy resolution data
    vector<double> sigPerMeV;   ///< gain calibration data
    double PEPerMeV;            ///< calibrated total PE per MeV
    bool ignoreOuter = false;   ///< ignore outer PMT ring in PMT sums
    
    /// set "ignore outer PMTs" mode
    void setIgnoreOuter(bool ig);
};

#endif
