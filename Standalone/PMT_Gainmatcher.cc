#include "SourceCalPlugin.hh"
#include "ReducedDataScanner.hh"
#include "OutputManager.hh"
#include "PathUtils.hh"
#include "StringManip.hh"
#include "HistogramSequenceFitter.hh"

#include <TStyle.h>

/// Interpolated backgrounds for AnalyzerPlugins
class PluginInterpolator {
public:
    /// Constructor
    PluginInterpolator() { }
    
    /// fit input data
    void fit();
    
    /// generate interpolated value
    void interpolate(RunAccumulator* RA) const;
    
    vector<RunAccumulator*> rundat;             ///< input data list
    
protected:
    ExponentialIntegralFitter IIF;              ///< core fitter routine
    vector<string> hNames;                      ///< names for each histogram
    vector<HistogramSequenceFitter> hFitters;   ///< fitters for each plugin histogram
    vector<intervalList> runIntervals;          ///< timing info for input data
    
    static intervalList getRunIntervals(const RunAccumulator* RA);
    static double t0;
};

double PluginInterpolator::t0 = 0;

intervalList PluginInterpolator::getRunIntervals(const RunAccumulator* RA) {
    intervalList L;
    if(!RA) return L;
    //double t0 = 0;
    for(auto it = RA->runTimes.counts.begin(); it != RA->runTimes.counts.end(); it++) {
        L.push_back(pair<double,double>(t0, t0 + it->second));
        t0 += it->second;
    }
    return L;
}

void PluginInterpolator::fit() {
    hNames.clear();
    hFitters.clear();
    runIntervals.clear();
    if(!rundat.size()) return;
    
    // collect histogram names
    const map<string,TH1*>& hlist = rundat[0]->getHists();
    for(auto it = hlist.begin(); it != hlist.end(); it++) hNames.push_back(it->first);
    // collect run intervals
    for(auto rit = rundat.begin(); rit != rundat.end(); rit++)
        runIntervals.push_back(getRunIntervals(*rit));
    // fit each histogram
    for(auto hnit = hNames.begin(); hnit != hNames.end(); hnit++) {
        hFitters.push_back(HistogramSequenceFitter(&IIF));
        for(size_t i = 0; i < rundat.size(); i++)
            hFitters.back().addData(rundat[i]->getSavedHist(*hnit),runIntervals[i]);
        hFitters.back().fit();
    }
}

void PluginInterpolator::interpolate(RunAccumulator* RA) const {
    if(!RA) return;
    intervalList L = getRunIntervals(RA);
    for(size_t i=0; i<hNames.size(); i++) {
        TH1* h = RA->getSavedHist(hNames[i]);
        hFitters[i].interpolate(L,h);
    }
}



class SourceRunSubtracter: public OutputManager, public PluginInterpolator {
public:
    SourceRunSubtracter(): OutputManager("Source_Calibrations",getEnvSafe("ACORN_SUMMARY")) {
        gStyle->SetOptStat("");
        defaultCanvas->SetLogy(true);
        defaultCanvas->SetLeftMargin(0.14);
    }
    
    void addBackgroundSegment(const vector<RunID>& rns) {
        SourceCalAnalyzer* SCAbg = new SourceCalAnalyzer(this);
        ReducedDataScanner Rb(false);
        Rb.addRuns(rns);
        SCAbg->loadProcessedData(Rb);
        rundat.push_back(SCAbg);
    }
    
    void analyzeForeground(const vector<RunID>& rns, const string& snm) {
        printf("Loading foreground data...\n");
        ReducedDataScanner Rf(false);
        Rf.addRuns(rns);
        SourceCalAnalyzer SCAfg(this);
        SCAfg.loadProcessedData(Rf);
        SCAfg.name += "_"+snm;
        SCAfg.mySourceCalPlugin->srcName = snm;
        
        printf("Interpolating background region...\n");
        SourceCalAnalyzer* SCAInterp = dynamic_cast<SourceCalAnalyzer*>(SCAfg.makeAnalyzer("interp_bg_"+snm, ""));
        SCAInterp->addSegment(SCAfg);
        interpolate(SCAInterp);
        
        SCAfg.mySourceCalPlugin->bgSubtrPlots(*SCAInterp->mySourceCalPlugin);
        SCAfg.setWriteRoot(true);
    }

};


void bg_subtraction_study() {
    
    SourceRunSubtracter SRS;
    
    printf("Generating background estimate...\n");

    for(unsigned int i=0; i<3; i++) {
        vector<RunID> bgruns;
        for(int j=0; j<4; j++) bgruns.push_back(RunID(1326+2*i,j));
        SRS.addBackgroundSegment(bgruns);
    }
    SRS.fit();
    
    vector<RunID> fgruns;
    for(int j=0; j<12; j++) fgruns.push_back(RunID(1327,j));
    SRS.analyzeForeground(fgruns,"Bi207");
    
    fgruns.clear();
    for(int j=0; j<10; j++) fgruns.push_back(RunID(1329,j));
    SRS.analyzeForeground(fgruns,"Sn113");
}

int main(int, char**) {
    bg_subtraction_study();
    return 0;
}
