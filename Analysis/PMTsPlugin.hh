#ifndef PMTSPLUGIN_HH
#define PMTSPLUGIN_HH

#include "RunAccumulator.hh"

/// Analyzer plugin for PMT (non-proton-coincident) data
class PMTsPlugin: public AnalyzerPlugin {
public:
    /// Constructor
    PMTsPlugin(RunAccumulator* RA);
    
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

/// Analyzer with wishbone plugin
class PMTsAnalyzer: public RunAccumulator {
public:
    /// Constructor
    PMTsAnalyzer(OutputManager* pnt, const std::string& nm = "PMTs", const std::string& inflName = ""): RunAccumulator(pnt,nm,inflName) {
        addPlugin(myPMTsPlugin = new PMTsPlugin(this));
    }
    
    /// create a new instance of this object (cloning self settings) for given directory
    virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) { return new PMTsAnalyzer(this,nm,inflname); }
    
    PMTsPlugin* myPMTsPlugin;   ///< PMTs analyzer plugin
};

#endif
