#include "SourceCalPlugin.hh"
#include "ReducedDataScanner.hh"
#include "OutputManager.hh"
#include "PathUtils.hh"
#include "strutils.hh"
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

void bg_subtraction_study() {
    
    OutputManager OM("Bi_PMT_Gain",getEnvSafe("ACORN_SUMMARY"));
    gStyle->SetOptStat("");
    OM.defaultCanvas->SetLogy(true);
    OM.defaultCanvas->SetLeftMargin(0.14);
    
    printf("Generating background estimate...\n");
    PluginInterpolator PI;
    SourceCalAnalyzer* SCAbg[3];
    for(unsigned int i=0; i<3; i++) {
        SCAbg[i] = new SourceCalAnalyzer(&OM);
        ReducedDataScanner Rb(false);
        for(int j=0; j<4; j++) Rb.addRun(RunID(1326+2*i,j));
        SCAbg[i]->loadProcessedData(Rb);
        PI.rundat.push_back(SCAbg[i]);
    }
    PI.fit();
       
    printf("Loading foreground data...\n");
    ReducedDataScanner Rf(false);
    for(int j=0; j<12; j++) Rf.addRun(RunID(1327,j));
    SourceCalAnalyzer SCAfg(&OM);
    SCAfg.loadProcessedData(Rf);
    
    printf("Interpolating background region...\n");
    SourceCalAnalyzer* SCAInterp = dynamic_cast<SourceCalAnalyzer*>(SCAfg.makeAnalyzer("interp_bg",""));
    SCAInterp->addSegment(SCAfg);
    PI.interpolate(SCAInterp);
    
    SCAfg.mySourceCalPlugin->bgSubtrPlots(*SCAInterp->mySourceCalPlugin);
    SCAfg.setWriteRoot(true);
}

int main(int, char**) {
    bg_subtraction_study();
    return 0;
}
