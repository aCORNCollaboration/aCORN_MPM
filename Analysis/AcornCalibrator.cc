#include "AcornCalibrator.hh"

AcornCalibrator::AcornCalibrator(RunID r): rn(r), md(MetadataDB::MDB.getRun(r)) {
}

double AcornCalibrator::calEnergy(double ADC) const {
    return md.calibPMT(ADC);
}
