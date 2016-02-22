/// \file RunAccumulator.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef RUNACCUMULATOR_HH
#define RUNACCUMULATOR_HH

#include "PluginSaver.hh"
#include "Enums.hh"
#include "TagCounter.hh"
#include "BaseDataScanner.hh"
#include "AcornDB.hh"
#include "SMFile.hh"

using std::vector;
using std::pair;

class RunAccumulatorPlugin;

enum DataMode {
    BAD,
    NG6,
    NGC
};

class RunAccumulator: public PluginSaver {
public:
    /// constructor
    RunAccumulator(OutputManager* pnt, const std::string& nm = "RunAccumulator", const std::string& inflName = "");
        
    /// add histograms from another RunAccumulator of the same type
    virtual void addSegment(const SegmentSaver& S) override;
    /// zero out run times
    void zeroCounters();
    /// scale all saved histograms by a factor
    virtual void scaleData(double s) override;
        
    /// write to SMFile
    virtual void write(string outName = "");
        
    /// fill data from a BaseDataScanner
    virtual void loadProcessedData(BaseDataScanner& BDS);
       
    bool isSimulated;               ///< flag for whether this is based on simulated data
    TagCounter<RunID> runCounts;    ///< event counts by run
    TagCounter<RunID> runTimes;     ///< time spent on each run
    SMFile qOut;
    
    /// fill core histograms in plugins from data point
    virtual void fillCoreHists(BaseDataScanner& PDS, double weight);
    /// calculate, upload analysis results
    virtual void makeAnaResults();
    /// run calculations and plots, save output files
    virtual void makeOutput(bool doPlots = true);
        
    /// merge every subdirectory of d containing analyzed data
    unsigned int mergeDir(const string& d);
    
    /// make rate-scaled copy of histogram; optionally divide by bin size on number of scale axes
    TH1* hToRate(TH1* h, int scaleAxes);
    /// generate base analysis result pre-filled with run range
    AnaResult makeBaseResult() const;
    
    DataMode dataMode = BAD;
    
protected:
    /// build plugins appropriate for input file; call in subclass after setting up myBuilders
    virtual void buildPlugins() override;
    vector<RunAccumulatorPlugin*> myRAPs;       ///< properly typecast active plugins
};

/// generic analyzer plug-in class for a RunAccumulator
class RunAccumulatorPlugin: public SegmentSaver {
public:
    /// constructor
    RunAccumulatorPlugin(RunAccumulator* RA, const string& nm, const string& inflname = ""):
    SegmentSaver(RA, nm, inflname), myA(RA) { }

    /// virtual routine for filling core histograms from data point
    virtual void fillCoreHists(BaseDataScanner& PDS, double weight) = 0;
    /// calculate, upload analysis results
    virtual void makeAnaResults() { }
    
    RunAccumulator* myA;    ///< RunAccumulator with which this plugin is associated
    /// shortcut for histogram to rate normalization
    TH1* hToRate(TH1* h, int scaleAxes) { return myA->hToRate(h, scaleAxes); }
};

/// Builder for RunAccumulatorPlugins
class RunAccumulatorPluginBuilder: public PluginBuilder {
public:
    /// Constructor
    RunAccumulatorPluginBuilder() { }
    
    virtual void makePlugin(OutputManager* pnt, const string& inflName = "") {
        RunAccumulator* RA = dynamic_cast<RunAccumulator*>(pnt);
        assert(RA);
        thePlugin = _makePlugin(RA,inflName);
    }
    
    virtual SegmentSaver* _makePlugin(RunAccumulator* RA, const string& inflName) = 0;
};

//////////////////////////////
//////////////////////////////

/// Helper class for range-based cuts
class RangeCutSet {
public:
    /// Constructor
    RangeCutSet() { }
    /// Start new region definition
    void newRegion() { regions.push_back(vector< pair<double,double> >()); norms.push_back(0); }
    /// Append range to newest region
    void appendRange(double x0, double x1) { 
        if(!regions.size()) newRegion();
        regions.back().push_back(pair<double,double>(x0,x1));
        norms.back() += x1-x0;
    }
    /// Scale all normalizations by given factor
    void scaleNorms(double s) { for(auto& n: norms) n *= s; }
    /// Scale all normalizations to specified
    void scaleNormsTo(size_t i) { scaleNorms(1./norms[i]); }
    
    vector< vector< pair<double,double> > > regions;    ///< the regions
    vector<double> norms;                               ///< normalization for each region
};

/// Helper class for background-region-subtracted histograms
class FGBGRegionsHist {
public:
    /// Constructor
    FGBGRegionsHist(RunAccumulatorPlugin* P): myP(P) { }
    /// Destructor
    ~FGBGRegionsHist() { for(auto h: hRates) delete h; }
    
    /// Setup using template histogram and regions
    void setTemplate(const TH1& hTemplate, const RangeCutSet& rs);    
   
    /// fill appropriate histogram
    void fill(double cutval, double x, double w);
    /// fill appropriate histogram, TH2 version
    void fill(double cutval, double x, double y, double w);
    
    /// generate rate-scaled, background-subtracted copies
    void makeRates(int axesScale = 1, double xscale = 1.0, bool bgsub = true);
    /// get regions list
    const RangeCutSet& getRegions() const { return myRegions; }
    /// mutable regions list access
    RangeCutSet& getRegions() { return myRegions; }
    
    vector<TH1*> hRates;        ///< background-subtracted, rate-scaled versions
    
protected:
    RunAccumulatorPlugin* myP;  ///< plugin to which this pair belongs
    vector<TH1*> hs;            ///< rate histograms for each region
    RangeCutSet myRegions;      ///< fill region definitions
};

#endif
