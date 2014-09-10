#ifndef WISHBONEPLUGIN_HH
#define WISHBONEPLUGIN_HH

#include "RunAccumulator.hh"
#include "ManualInfo.hh"

#include <vector>
#include <utility>
using std::vector;
using std::pair;

/// Utility class for wishbone plot backrgound subtraction and projections
class TH2Slicer {
public:
    /// Constructor
    TH2Slicer(TH2* hh): h(hh) { }
    
    /// Project contents in bins closest to yrange; return actual interval summed
    double projSlice(double y0, double y1, TH1*& projOut);
    
    /// Make a background-profile-subtracted copy of a histogram
    TH2* subtractProfile(const TH1* p, double s=1.) const;
        
protected:
    static int pcount;        ///< projection number counter for unique naming
    TH2* h;
};

/// Analyzer plugin for wishbone data
class WishbonePlugin: public AnalyzerPlugin {
public:
    /// Constructor
    WishbonePlugin(RunAccumulator* RA);
    
    /// Fill core histograms from data point
    virtual void fillCoreHists(BaseDataScanner& PDS, double weight);
    
    /// generate calculated hists
    virtual void calculateResults();
    /// Generate output plots
    virtual void makePlots();
    
    
    TH1* hProtonSignal[2];      ///< proton detector signal, electron coincident or not
    TH1* h4pTiming;             ///< "4p" mode discriminator arrival time
    TH2F* hWishbone;            ///< wishbone plot
    TH1* hNVeto;                ///< number of veto PMTs for wishbone-like events
    TH2F* hNE;                  ///< number of main PMTs triggered vs. total event energy
    TH2F* hPMTs;                ///< which PMTs fired, as a function of event
    TH2F* hChanSpec;            ///< individual PMTs spectrum distribution
    TH2F* hModuleMult;          ///< module multiplicity, for module dropout issues
    
    double E_p_lo = 650;        ///< proton signal low cut for wishbone data
    double E_p_hi = 2200;       ///< proton signal high cut for wishbone data
    double T_p_min = 750;       ///< minimum TOF for background analysis [ns]
    double T_p_lo = 2750;       ///< proton TOF lower window for wishbone [ns]
    double T_p_hi = 4500;       ///< proton TOF lower window for wishbone [ns]
    double T_p_max = 9500;      ///< maximum TOF for background analysis [ns]
    
    TH1* hWishboneEProj[2];     ///< Wishbone energy spectrum, background and background-subtracted
    TH1* hWishboneTProj;        ///< Wishbone time-axis projection
    TH2* hWishboneBGSub;        ///< Background-subtracted wishbone
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
