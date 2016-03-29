/// \file PMT_Gainmatcher.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "SourceCalPlugin.hh"
#include "PluginInterpolator.hh"
#include "ReducedDataScanner.hh"
#include "PathUtils.hh"
#include "StringManip.hh"

#include <TStyle.h>
#include <cassert>

class SourceCalAnalyzer: public RunAccumulator {
public:
    SourceCalAnalyzer(OutputManager* pnt, const std::string& nm = "SourceCal", const std::string& inflname = ""):
    RunAccumulator(pnt, nm, inflname) {
        myBuilders["SourceCalPlugin"] = &mySourceCalPluginBuilder;
        buildPlugins();
    }
    
    RunAccumulatorPluginBuilder<SourceCalPlugin> mySourceCalPluginBuilder;
};

class SourceRunSubtracter: public OutputManager, public PluginInterpolator {
public:
    SourceRunSubtracter(): OutputManager("Source_Calibrations",getEnvSafe("ACORN_SUMMARY")) {
        gStyle->SetOptStat("");
        defaultCanvas.SetLogy(true);
        defaultCanvas.SetLeftMargin(0.14);
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
        dynamic_cast<SourceCalPlugin*>(SCAfg.mySourceCalPluginBuilder.thePlugin)->srcName = snm;
        
        printf("Interpolating background region...\n");
        SourceCalAnalyzer* SCAInterp = dynamic_cast<SourceCalAnalyzer*>(SCAfg.makeAnalyzer("interp_bg_"+snm, ""));
        SCAInterp->addSegment(SCAfg);
        interpolate(SCAInterp);
        
        //dynamic_cast<SourceCalPlugin*>(SCAfg.mySourceCalPluginBuilder.thePlugin)->bgSubtrPlots(*dynamic_cast<SourceCalPlugin*>(SCAInterp.mySourceCalPluginBuilder.thePlugin));
        assert(false); // TODO line above
        SCAfg.writeROOT();
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
    
    // Bi207 data
    vector<RunID> fgruns;
    for(int j=0; j<12; j++) fgruns.push_back(RunID(1327,j));
    SRS.analyzeForeground(fgruns,"Bi207");
    
    // Sn113 data
    fgruns.clear();
    for(int j=0; j<10; j++) fgruns.push_back(RunID(1329,j));
    SRS.analyzeForeground(fgruns,"Sn113");
}

int main(int, char**) {
    bg_subtraction_study();
    return 0;
}
