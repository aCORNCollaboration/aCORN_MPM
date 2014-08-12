#ifndef WISHBONEPLUGIN_HH
#define WISHBONEPLUGIN_HH

#include "RunAccumulator.hh"
#include "ManualInfo.hh"

/// Analyzer plugin for wishbone data
class WishbonePlugin: public AnalyzerPlugin {
public:
    /// Constructor
    WishbonePlugin(RunAccumulator* RA);
    
    /// Fill core histograms from data point
    virtual void fillCoreHists(BaseDataScanner& PDS, double weight);
    
    /// Generate output plots
    virtual void makePlots();
    
    TH1* hProtonSignal;         ///< proton detector signal
    TH2F* hWishbone;            ///< wishbone plot
};


/// Analyzer with wishbone plugin
class WishboneAnalyzer: public RunAccumulator {
public:
    /// Constructor
    WishboneAnalyzer(OutputManager* pnt, const std::string& nm = "Wishbone", const std::string& inflName = ""): RunAccumulator(pnt,nm,inflName) {
        addPlugin(myWishbonePlugin = new WishbonePlugin(this));
    }
    
    /// create a new instance of this object (cloning self settings) for given directory
    virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) { return new WishboneAnalyzer(this,nm,inflname); }
    
    WishbonePlugin* myWishbonePlugin;   ///< Wishbone analyzer plugin
};

#endif
