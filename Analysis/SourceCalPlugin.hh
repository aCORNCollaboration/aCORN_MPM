/// \file SourceCalPlugin.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef SOURCECALPLUGIN
#define SOURCECALPLUGIN

#include "RunAccumulator.hh"

/// Analyzer plugin for source calibrations (non-proton-coincident data)
class SourceCalPlugin: public RunAccumulatorPlugin {
public:
    /// Constructor
    SourceCalPlugin(RunAccumulator* RA, const string& nm = "SourceCalPlugin", const string& inflname = "");
    
    /// Fill core histograms from data point
    void fillCoreHists(BaseDataScanner& PDS, double weight) override;
    
    /// Generate background-subtracted output
    void bgSubtrPlots(SourceCalPlugin& bg);
    
    TH1* hEnergy;               ///< Calibrated energy spectrum
    TH1* hEnergyRecal;          ///< Individual-PMT-based recalibrated spectrum
    TH1* hPMTSig[N_E_PMT];      ///< raw PMT signal spectra
    
    string srcName;             ///< type of source being analyzed
};

#endif
