/// \file AcornCalibrator.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "AcornCalibrator.hh"
#include "AcornDB.hh"
#include "BaseDataScanner.hh"

AcornCalibrator::AcornCalibrator(RunID r): rn(r) {
    printf("Initializing energy calibrator for %i:%i",r.first,r.second);
    
    AcornDB::ADB().getPMTSumCal(rn, slope, intercept);
    AcornDB::ADB().getPMTcal(rn, sigPerPE, sigPerMeV);
    gRecal = AcornDB::ADB().getRecal(rn);
    
    setIgnoreOuter(false);
}

AcornCalibrator::~AcornCalibrator() {
    if(gRecal) delete gRecal;
}

void AcornCalibrator::setIgnoreOuter(bool ig) {
    ignoreOuter = ig;
    
    PEPerMeV = 0;
    for(size_t i=0; i<N_E_PMT; i++) {
        if(ignoreOuter && isOuter(i)) continue;
        PEPerMeV += sigPerMeV[i] / sigPerPE[i];
    }
    printf(" (%.1f PE/MeV; %g + %g*s%s).\n",PEPerMeV,intercept,slope,gRecal?" +Recal":"");
}

bool AcornCalibrator::isOuter(size_t i) {
    if(i<=3 || i >= 15) return true;
    return i==6 || i==7 || i==11 || i==12;
}

double AcornCalibrator::wsum(const Short_t* ADC) const {
    double npe = 0;
    for(size_t i=0; i<N_E_PMT; i++) {
        if(ignoreOuter && isOuter(i)) continue;
        npe += ADC[i] / sigPerPE[i];
    }
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

double AcornCalibrator::calWishboneSum(double ADC) const {
    if(!gRecal) return calPMTSum(ADC);
    return gRecal->Eval(slope*ADC) + intercept;
}
