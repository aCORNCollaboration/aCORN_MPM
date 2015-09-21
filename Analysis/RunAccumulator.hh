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

class RunAccumulator: public PluginSaver {
public:
    /// constructor
    RunAccumulator(OutputManager* pnt, const std::string& nm = "RunAccumulator", const std::string& inflName = "");
        
    /// add histograms from another RunAccumulator of the same type
    virtual void addSegment(const SegmentSaver& S);
    /// zero out run times
    void zeroCounters();
    /// scale all saved histograms by a factor
    virtual void scaleData(double s);
        
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
        
    /// create a new instance of this object (cloning self settings) for given directory
    virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) { return new RunAccumulator(this,nm,inflname); }
        
    /// merge every subdirectory of basePath containing analyzed data
    unsigned int mergeDir();
    
    /// make rate-scaled copy of histogram; optionally divide by bin size on number of scale axes
    TH1* hToRate(TH1* h, int scaleAxes);
    /// generate base analysis result pre-filled with run range
    AnaResult makeBaseResult() const;
    
protected:
    /// build plugins appropriate for input file; call in subclass after setting up myBuilders
    virtual void buildPlugins();
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

/// Helper class for background-region-subtracted histograms
class FGBGRegionsHist {
public:
    /// Constructor
    FGBGRegionsHist(RunAccumulatorPlugin* P);
    /// Destructor
    ~FGBGRegionsHist();
    /// Setup using template histogram
    void setTemplate(const TH1& hTemplate);    
    /// add region to count as FG/BG
    void addRegion(double x0, double x1, bool fg);
    /// fill appropriate histogram
    void fill(double cutval, double x, double w);
    /// fill appropriate histogram, TH2 version
    void fill(double cutval, double x, double y, double w);
    /// generate rate-scaled, background-subtracted copies
    void makeRates(int axesScale = 1, double xscale = 1.0);

    TH1* hRates[2];     ///< background-subtracted, rate-scaled versions
    
protected:
    RunAccumulatorPlugin* myP;                  ///< plugin to which this pair belongs
    TH1* h[2];                                  ///< background, foreground region histograms
    vector< pair<double,double> > regions[2];   ///< regions defined for FG/BG
    double totalLength[2];                      ///< total length of FG/BG regions
};

#endif