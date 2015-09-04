/// \file PMTsPlugin.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef PMTSPLUGIN_HH
#define PMTSPLUGIN_HH

#include "RunAccumulator.hh"

/// Analyzer plugin for PMT (non-proton-coincident) data
class PMTsPlugin: public RunAccumulatorPlugin {
public:
    /// Constructor
    PMTsPlugin(RunAccumulator* RA, OutputManager* pnt, const string& nm, const string& inflname = "");
    
    /// Fill core histograms from data point
    virtual void fillCoreHists(BaseDataScanner& PDS, double weight);
    
    /// generate calculated hists
    //virtual void calculateResults();
    /// Generate output plots
    virtual void makePlots();
   
    TH1* hEnergy;               ///< event energy spectrum
    TH2* hChanSpec;             ///< individual PMTs spectrum distribution
    TH2* hNE;                   ///< number of main PMTs triggered vs. total event energy
    TH1* hEETime;               ///< timing between electron events
    
protected:
    
    double prev_e_time;         ///< time of previous electron event
};

/// Builder for RunAccumulatorPlugins
class PMTsPluginBuilder: public RunAccumulatorPluginBuilder {
public:
    /// Constructor
    PMTsPluginBuilder(RunAccumulator* R): RunAccumulatorPluginBuilder(R) { }
    /// instantiate plugin SegmentSaver
    virtual void makePlugin(OutputManager* pnt, const string& inflName = "") { thePlugin = new PMTsPlugin(RA, pnt, "PMTs", inflName); }
};

#endif
