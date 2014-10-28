#include "AcornCalibrator.hh"
#include "AcornDB.hh"

AcornCalibrator::AcornCalibrator(RunID r): rn(r) {
    AcornDB::ADB().getPMTSumCal(rn, slope, intercept);
}

double AcornCalibrator::calEnergy(double ADC) const {
    return slope*ADC + intercept;
}
