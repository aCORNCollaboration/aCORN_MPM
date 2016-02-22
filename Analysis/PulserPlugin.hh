/// \file PulserPlugin.hh aCORN proton pulser signal analysis plugin
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2016

#ifndef PULSERPLUGIN_HH
#define PULSERPLUGIN_HH

#include "RunAccumulator.hh"
#include <deque>
using std::deque;

/// Analyzer plugin for periodic pulser data
class PulserPlugin: public RunAccumulatorPlugin {
public:
    /// Constructor
    PulserPlugin(RunAccumulator* RA, const string& nm, const string& inflname = "");
    
    /// Fill core histograms from data point
    void fillCoreHists(BaseDataScanner& PDS, double weight) override;
    /// generate calculated hists
    void calculateResults() override;
    /// calculate, upload analysis results
    void makeAnaResults() override;
    /// Generate output plots
    void makePlots() override;
    
    double E_puls_lo = 7000;            ///< proton pulser low cut
    double E_puls_hi = 12000;           ///< proton pulser high cut
    double T_pulser = 0.9958e9;         ///< expected pulser timing [ns]
    double tau_long = 0.05;             ///< threshold for identifying dropouts [ns]
    
    deque<double> prevPulses;           ///< times of preceding pulses
    double T_prev_p = 0;                ///< time of preceding proton event [ns]
    double E_prev_p = 0;                ///< detector signal of preceding proton [ADC channels]
    FGBGRegionsHist hPulserSignal;      ///< proton detector signal, cut on pulser timing
    TH1* hPulserTiming;                 ///< periodic pulser timing, broad window
    TH1* hPulserTimingPrecis;           ///< periodic pulser timing, precision window
    TH1* hPTimingGap;                   ///< time gap between subsequent proton events
    TVectorD* vGapTime;                 ///< DAQ deadtime inferred from large "gaps" [s]
    FGBGRegionsHist hTimingHumpPrev;    ///< proton detector signal in correlated events hump versus "normal" peak
    FGBGRegionsHist hTimingHumpPost;    ///< proton detector signal in correlated events hump versus "normal" peak
    
protected:
    /// set up analysis cuts for data mode
    void setAnalysisCuts();
    RangeCutSet psignalRegions;         ///< proton signal energy regions
    RangeCutSet dtpRegions;             ///< proton delta time regions
};

/// Builder for PulserPlugins
class PulserPluginBuilder: public RunAccumulatorPluginBuilder {
public:
    /// Constructor
    PulserPluginBuilder() { }
    /// instantiate plugin SegmentSaver
    virtual SegmentSaver* _makePlugin(RunAccumulator* RA, const string& inflName) override { return new PulserPlugin(RA, "PulserPlugin", inflName); }
};

#endif
