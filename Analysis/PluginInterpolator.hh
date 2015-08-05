/// \file PluginInterpolator.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef PLUGININTERPOLATOR_HH
#define PLUGININTERPOLATOR_HH

#include "RunAccumulator.hh"
#include "HistogramSequenceFitter.hh"

/// Interpolated backgrounds for AnalyzerPlugins
class PluginInterpolator {
public:
    /// Constructor
    PluginInterpolator() { }
    
    /// fit input data
    void fit();
    
    /// generate interpolated value
    void interpolate(RunAccumulator* RA) const;
    
    vector<RunAccumulator*> rundat;             ///< input data list
    
protected:
    ExponentialIntegralFitter IIF;              ///< core fitter routine
    vector<string> hNames;                      ///< names for each histogram
    vector<HistogramSequenceFitter> hFitters;   ///< fitters for each plugin histogram
    vector<intervalList> runIntervals;          ///< timing info for input data
    
    static intervalList getRunIntervals(const RunAccumulator* RA);
    static double t0;
};

#endif
