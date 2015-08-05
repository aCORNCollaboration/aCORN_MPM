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
class SourceCalPlugin: public AnalyzerPlugin {
public:
    /// Constructor
    SourceCalPlugin(RunAccumulator* RA);
    
    /// Fill core histograms from data point
    virtual void fillCoreHists(BaseDataScanner& PDS, double weight);
    
    /// generate calculated hists
    virtual void calculateResults() { }
    /// Generate output plots
    virtual void makePlots() { }
    /// Generate background-subtracted output
    void bgSubtrPlots(SourceCalPlugin& bg);
    
    TH1* hEnergy;               ///< Calibrated energy spectrum
    TH1* hEnergyRecal;          ///< Individual-PMT-based recalibrated spectrum
    TH1* hPMTSig[N_E_PMT];      ///< raw PMT signal spectra
    
    string srcName;             ///< type of source being analyzed
};


/// Analyzer with SourceCal plugin
class SourceCalAnalyzer: public RunAccumulator {
public:
    /// Constructor
    SourceCalAnalyzer(OutputManager* pnt, const std::string& nm = "SourceCal", const std::string& inflName = ""): RunAccumulator(pnt,nm,inflName) {
        addPlugin(mySourceCalPlugin = new SourceCalPlugin(this));
    }
    
    /// create a new instance of this object (cloning self settings) for given directory
    virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) { return new SourceCalAnalyzer(this,nm,inflname); }
    
    SourceCalPlugin* mySourceCalPlugin;         ///< Source calibration analyzer plugin
};


#endif
