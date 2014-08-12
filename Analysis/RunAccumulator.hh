#ifndef RUNACCUMULATOR_HH
#define RUNACCUMULATOR_HH

#include "SegmentSaver.hh"
#include "Enums.hh"
#include "TagCounter.hh"
#include "BaseDataScanner.hh"

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
        
        /// write to QFile
        virtual void write(std::string outName = "");
        
        /// fill data from a BaseDataScanner
        virtual void loadProcessedData(BaseDataScanner& BDS);
       
        bool isSimulated;               ///< flag for whether this is based on simulated data
        TagCounter<RunID> runCounts;    ///< event counts by run
        TagCounter<RunID> runTimes;     ///< time spent on each run
        
        /// fill core histograms in plugins from data point
        virtual void fillCoreHists(BaseDataScanner& PDS, double weight);
        /// calculate results from filled histograms
        virtual void calculateResults();
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
                
protected:

        std::map<std::string,AnalyzerPlugin*> myPlugins;        ///< analysis plugins
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
        void printCanvas(std::string fname, std::string suffix=".pdf") const { myA->printCanvas(fname,suffix); }
        
        std::string name;       ///< plugin name
        RunAccumulator* myA;    ///< RunAccumulator with which this plugin is associated
        
        /// virtual routine for filling core histograms from data point
        virtual void fillCoreHists(BaseDataScanner& PDS, double weight) = 0;
        
        /// generate output plots
        virtual void makePlots() {}
        /// generate calculated hists
        virtual void calculateResults() {}
        /// virtual routine for MC/Data comparison plots/calculations
        /// NOTE: this MUST NOT change the contents of saved histograms (calculated ones are OK)
        virtual void compareMCtoData(AnalyzerPlugin*) {}

protected:
        std::vector<TH1*> myHists; ///< histograms registered by this plugin
};

#endif