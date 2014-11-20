#include "AcornCalibrator.hh"
#include "AcornDB.hh"
#include "BaseDataScanner.hh"

AcornCalibrator::AcornCalibrator(RunID r): rn(r) {
    printf("Initializing energy calibrator for %i:%i",r.first,r.second);
    
    AcornDB::ADB().getPMTSumCal(rn, slope, intercept);
    AcornDB::ADB().getPMTcal(rn, sigPerPE, sigPerMeV);
    
    PEPerMeV = 0;
    for(size_t i=0; i<N_E_PMT; i++) PEPerMeV += sigPerMeV[i] / sigPerPE[i];
    printf(" (%.1f PE/MeV; %g + %g*s).\n",PEPerMeV,intercept,slope);
}

double AcornCalibrator::wsum(const Short_t* ADC) const {
    double npe = 0;
    for(size_t i=0; i<N_E_PMT; i++) npe += ADC[i] / sigPerPE[i];
    return npe;
}

void AcornCalibrator::invcalPMTSum(CalPeak& pk) const {
    double s =  invcalPMTSum(pk.center.x);
    double d = dEds(s);
    pk.center.x = s;
    pk.center.err /= d;
    pk.sigma = (1/d)*pk.sigma;
    pk.height = d*pk.height;
}

void AcornCalibrator::invcalWSum(CalPeak& pk) const {
    double w =  invcalWSum(pk.center.x);
    double d = dEdw(w);
    pk.center.x = w;
    pk.center.err /= d;
    pk.sigma = (1/d)*pk.sigma;
    pk.height = d*pk.height;
}
