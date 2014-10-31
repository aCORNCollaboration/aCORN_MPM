#include "AcornCalibrator.hh"
#include "AcornDB.hh"
#include "BaseDataScanner.hh"

AcornCalibrator::AcornCalibrator(RunID r): rn(r) {
    printf("Initializing energy calibrator for %i:%i",r.first,r.second);
    
    AcornDB::ADB().getPMTSumCal(rn, slope, intercept);
    AcornDB::ADB().getPMTcal(rn, sigPerPE, sigPerMeV);
    
    PEPerMeV = 0;
    for(size_t i=0; i<N_E_PMT; i++) PEPerMeV += sigPerMeV[i] / sigPerPE[i];
    printf(" (%.1f PE/MeV).\n",PEPerMeV);
}

double AcornCalibrator::calEnergy(double ADC) const {
    return slope*ADC + intercept;
}

double AcornCalibrator::calEnergy(const Short_t* ADC) const {
    double npe = 0;
    for(size_t i=0; i<N_E_PMT; i++) npe += ADC[i] / sigPerPE[i];
    return 1000 * npe / PEPerMeV;
}
