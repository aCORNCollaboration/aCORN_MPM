#ifndef RUNACCUMULATOR_HH
#define RUNACCUMULATOR_HH

#include "SegmentSaver.hh"
#include "Enums.hh"
#include "TagCounter.hh"
#include "BaseDataScanner.hh"
#include "AcornDB.hh"

using std::vector;
using std::pair;

class AnalyzerPlugin;

class RunAccumulator: public SegmentSaver {
public:
    /// constructor
    RunAccumulator(OutputManager* pnt, const std::string& nm = "RunAccumulator", const std::string& inflName = "");
    /// destructor
    virtual ~RunAccumulator();
        
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
        
    /// fill core histograms in plugins from data point
    virtual void fillCoreHists(BaseDataScanner& PDS, double weight);
    /// calculate results from filled histograms
    virtual void calculateResults();
    /// calculate, upload analysis results
    virtual void makeAnaResults();
    /// make plots from each plugin
    virtual void makePlots();
    /// run calculations and plots, save output files
    virtual void makeOutput(bool doPlots = true);
    /// MC/Data comparison plots/calculations from each plugin
    virtual void compareMCtoData(RunAccumulator& OAdata);
    /// add an analyzer plugin
    void addPlugin(AnalyzerPlugin* AP);
    /// get plugin by name
    AnalyzerPlugin* getPlugin(const std::string& nm);
        
    /// create a new instance of this object (cloning self settings) for given directory
    virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) { return new RunAccumulator(this,nm,inflname); }
        
    /// merge every subdirectory of basePath containing analyzed data
    unsigned int mergeDir();
    
    /// make rate-scaled copy of histogram; optionally divide by bin size on number of scale axes
    TH1* hToRate(TH1* h, int scaleAxes);
    /// generate base analysis result pre-filled with run range
    AnaResult makeBaseResult() const;
    
protected:
    map<std::string,AnalyzerPlugin*> myPlugins;        ///< analysis plugins
};

/// generic analyzer plug-in class
class AnalyzerPlugin {
public:
        /// constructor
        AnalyzerPlugin(RunAccumulator* RA, const std::string& nm): name(nm), myA(RA) { }
        /// destructor
        virtual ~AnalyzerPlugin() {}
        
        /// create or load a TH1F stored histogram
        TH1* registerHist(const std::string& hname, const std::string& title, unsigned int nbins, float xmin, float xmax);
        /// create or load a stored histogram from template
        TH1* registerHist(const std::string& nm, const TH1& hTemplate);
        /// save canvas image
        void printCanvas(string fname, string suffix=".pdf") const { myA->printCanvas(fname,suffix); }
        
        string name;            ///< plugin name
        RunAccumulator* myA;    ///< RunAccumulator with which this plugin is associated
        
        /// virtual routine for filling core histograms from data point
        virtual void fillCoreHists(BaseDataScanner& PDS, double weight) = 0;
        
        /// generate calculated/normalized histograms derived from core saved data; possibly needed below.
        virtual void calculateResults() {}
        /// calculate, upload analysis results
        virtual void makeAnaResults() {}
        /// generate output plots
        virtual void makePlots() {}
        /// virtual routine for MC/Data comparison plots/calculations
        /// NOTE: this MUST NOT change the contents of saved histograms (calculated ones are OK)
        virtual void compareMCtoData(AnalyzerPlugin*) {}

protected:
        vector<TH1*> myHists; ///< histograms registered by this plugin
};

/// Helper class for background-region-subtracted histograms
class FGBGRegionsHist {
public:
    /// Constructor
    FGBGRegionsHist(AnalyzerPlugin* P);
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
    AnalyzerPlugin* myP;                        ///< plugin to which this pair belongs
    TH1* h[2];                                  ///< background, foreground region histograms
    vector< pair<double,double> > regions[2];   ///< regions defined for FG/BG
    double totalLength[2];                      ///< total length of FG/BG regions
};

#endif