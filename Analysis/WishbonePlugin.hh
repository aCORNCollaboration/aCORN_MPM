/// \file WishbonePlugin.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef WISHBONEPLUGIN_HH
#define WISHBONEPLUGIN_HH

#include "RunAccumulator.hh"
#include "TRateCategories.hh"

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
class WishbonePlugin: public RunAccumulatorPlugin {
public:
    /// Constructor
    WishbonePlugin(RunAccumulator* RA, const string& nm, const string& inflname = "");
    
    /// Fill core histograms from data point
    void fillCoreHists(BaseDataScanner& PDS, double weight) override;
    
    /// generate calculated hists
    void calculateResults() override;
    /// calculate, upload analysis results
    void makeAnaResults() override;
    /// Generate output plots
    void makePlots() override;
    
    
    FGBGRegionsHist hProtonSignal;      ///< proton detector signal
    //TH1* h4pTiming;                     ///< "4p" mode discriminator arrival time
    TH2F* hWishbone;                    ///< wishbone counts
    TH2F* hWBRate;                      ///< wishbone, non-background-subtracted, converted to rate
    FGBGRegionsHist hNVeto;             ///< number of veto PMTs for wishbone-like events
    FGBGRegionsHist hVetoSum;           ///< veto PMTs signal sum
    FGBGRegionsHist hNE;                ///< number of main PMTs triggered vs. total event energy
    FGBGRegionsHist hPMTs;              ///< which PMTs fired, as a function of event
    FGBGRegionsHist hChanSpec;          ///< individual PMTs spectrum distribution
    FGBGRegionsHist hModuleMult;        ///< module multiplicity, for module dropout issues
    
    /// types of tabulated rates
    enum RateType {
        RATE_WB_BG = 0,                 ///< Wishbone background electron/proton coincidences
        RATE_WB_FG = 1,                 ///< Wishbone foreground electron/proton coincidences
        RATE_PROTON = 2,                ///< All proton detector events
        RATE_VETO_BG = 3,               ///< veto-triggered background events
        RATE_VETO_FG = 4                ///< veto-triggered foreground events
    };
    TRateCategories* hRateHistory;      ///< Tabulated event rates
    
    FGBGRegionsHist hPos;               ///< PMT hit center-of-mass positions
    FGBGRegionsHist hPosSigma;          ///< PMT hit RMS spread
    FGBGRegionsHist hEnergyRadius;      ///< PMT hit position radius vs energy
    
    double E_p_lo = 650;                ///< proton signal low cut for wishbone data
    double E_p_hi = 2400;               ///< proton signal high cut for wishbone data
    double T_p_min = 750;               ///< minimum TOF for background analysis [ns]
    double T_p_lo = 2750;               ///< proton TOF lower window for wishbone [ns]
    double T_p_hi = 4500;               ///< proton TOF lower window for wishbone [ns]
    double T_p_max = 9500;              ///< maximum TOF for background analysis [ns]
    double Escale = 1.0;                ///< calculated energy re-scaling factor
        
    TH1* hWishboneEProj[2];             ///< Wishbone energy spectrum, background and background-subtracted
    TH1* hWishboneTProj;                ///< Wishbone time-axis projection [Hz/us]
    TH2* hWishboneBGSub;                ///< Background-subtracted wishbone [Hz/MeV/us]
    TH1* hWishboneFiducialTProj;        ///< Background-subtracted time profile in 100--300 keV region
    
protected:
    /// set up analysis cuts for data mode
    void setAnalysisCuts();
    /// initialize foreground/background regions
    void initRegions(FGBGRegionsHist& h);
};

/// Builder for RunAccumulatorPlugins
class WishbonePluginBuilder: public RunAccumulatorPluginBuilder {
public:
    /// Constructor
    WishbonePluginBuilder() { }
    /// instantiate plugin SegmentSaver
    virtual SegmentSaver* _makePlugin(RunAccumulator* RA, const string& inflName) override { return new WishbonePlugin(RA, "WishbonePlugin", inflName); }
};

#include "PulserPlugin.hh"

/// Analyzer class using WishbonePlugin
class WishboneAnalyzer: public RunAccumulator {
public:
    WishboneAnalyzer(OutputManager* pnt, const std::string& nm = "Wishbone", const std::string& inflname = ""):
    RunAccumulator(pnt, nm, inflname) {
        myBuilders["WishbonePlugin"] = &myWishbonePluginBuilder;
        myBuilders["PulserPlugin"] = &myPulserPluginBuilder;
        buildPlugins();
    }
    
    /// create a new instance of this object (cloning self settings) for given directory
    virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) override { return new WishboneAnalyzer(this,nm,inflname); }
    
    WishbonePluginBuilder myWishbonePluginBuilder;
    PulserPluginBuilder myPulserPluginBuilder;
};

#endif
