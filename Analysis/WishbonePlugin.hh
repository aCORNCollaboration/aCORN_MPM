/// \file WishbonePlugin.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef WISHBONEPLUGIN_HH
#define WISHBONEPLUGIN_HH

#include "RunAccumulator.hh"
#include "RangeConfigFile.hh"
#include "AcornDB.hh"

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
class WishbonePlugin: public RunAccumulatorPlugin {
public:
    /// Constructor
    WishbonePlugin(RunAccumulator* RA, OutputManager* pnt, const string& nm, const string& inflname = "");
    
    /// Fill core histograms from data point
    virtual void fillCoreHists(BaseDataScanner& PDS, double weight);
    
    /// generate calculated hists
    virtual void calculateResults();
    /// calculate, upload analysis results
    virtual void makeAnaResults();
    /// Generate output plots
    virtual void makePlots();
    
    
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
    
    FGBGRegionsHist hPos;               ///< PMT hit center-of-mass positions
    FGBGRegionsHist hPosSigma;          ///< PMT hit RMS spread
    FGBGRegionsHist hEnergyRadius;      ///< PMT hit position radius vs energy
    
    double E_p_lo = 650;                ///< proton signal low cut for wishbone data
    double E_p_hi = 2400;               ///< proton signal high cut for wishbone data
    double T_p_min = 750;               ///< minimum TOF for background analysis [ns]
    double T_p_lo = 2750;               ///< proton TOF lower window for wishbone [ns]
    double T_p_hi = 4500;               ///< proton TOF lower window for wishbone [ns]
    double T_p_max = 9500;              ///< maximum TOF for background analysis [ns]
    
    TH1* hWishboneEProj[2];             ///< Wishbone energy spectrum, background and background-subtracted
    TH1* hWishboneTProj;                ///< Wishbone time-axis projection [Hz/us]
    TH2* hWishboneBGSub;                ///< Background-subtracted wishbone [Hz/MeV/us]
    TH1* hWishboneFiducialTProj;        ///< Background-subtracted time profile in 100--300 keV region

protected:
    /// initialize foreground/background regions
    void initRegions(FGBGRegionsHist& h);
    
    /// set NG-C cuts
    void config_NGC_cuts();
};

/// Builder for RunAccumulatorPlugins
class WishbonePluginBuilder: public RunAccumulatorPluginBuilder {
public:
    /// Constructor
    WishbonePluginBuilder(RunAccumulator* R): RunAccumulatorPluginBuilder(R) { }
    /// instantiate plugin SegmentSaver
    virtual void makePlugin(OutputManager* pnt, const string& inflName = "") { thePlugin = new WishbonePlugin(RA, pnt, "Wishbone", inflName); }
};

#endif
