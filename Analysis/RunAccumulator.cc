/// \file RunAccumulator.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "RunAccumulator.hh"
#include "GraphUtils.hh"
#include "SMExcept.hh"
#include "StringManip.hh"
#include "PathUtils.hh"
#include <time.h>
#include <TH1.h>

RunAccumulator::RunAccumulator(OutputManager* pnt, const std::string& nm, const std::string& inflName):
PluginSaver(pnt,nm,inflName), isSimulated(false) {
    dataMode = BAD;
    string dsrc = split(strip(getEnvSafe("ACORN_REDUCED_ROOT"),"/"),"/").back();
    if(dsrc == "ROOT_NG6") dataMode = NG6;
    if(dsrc == "ROOT_NGC") dataMode = NGC;
    assert(dataMode != BAD);
    
    TH1::SetDefaultSumw2(true); // all histograms default to having errorbars
        
    TCumulativeMap<RunNum,Double_t> TCMTemplate("TCMTemplate");
    runCounts = (TCumulativeMap<RunNum,Double_t>*)registerCumulative("runCounts", TCMTemplate);
    runTimes = (TCumulativeMap<RunNum,Double_t>*)registerCumulative("runTimes", TCMTemplate);
}

void RunAccumulator::buildPlugins() {
    PluginSaver::buildPlugins();
    /// identify RunAccumulatorPlugins
    for(auto& kv: myBuilders) {
        auto P = dynamic_cast<RunAccumulatorPlugin*>(kv.second->thePlugin.get());
        if(P) myRAPs.push_back(P);
    }
}

void RunAccumulator::fillCoreHists(BaseDataScanner& PDS, double weight) {
    for(auto RAP: myRAPs) RAP->fillCoreHists(PDS,weight);
}

AnaResult RunAccumulator::makeBaseResult() const {
    AnaResult baseResult;
    if(runCounts->Size()) {
        baseResult.start = toRunID(runCounts->GetData().begin()->first);
        baseResult.end = toRunID(runCounts->GetData().rbegin()->first);
    }
    return baseResult;
}

void RunAccumulator::makeAnaResults() {
    if(!isCalculated) calculateResults();
    AcornDB::ADB().beginTransaction();
    
    AnaResult baseResult = makeBaseResult();
    baseResult.value = runCounts->GetTotal();
    AcornDB::ADB().uploadAnaResult("total_counts", "Total analyzed counts", baseResult);
    baseResult.value = runTimes->GetTotal();
    AcornDB::ADB().uploadAnaResult("total_time", "Total analyzed time [s]", baseResult);
    for(auto RAP: myRAPs) RAP->makeAnaResults();
    AcornDB::ADB().endTransaction();
}

void RunAccumulator::loadProcessedData(BaseDataScanner& PDS) {
    printf("Loading processed data...\n");
    if(!PDS.getnFiles())
        return;
    PDS.startScan();
    PDS.T_0 = 0;
    unsigned int nScanned = 0;
    
    while(PDS.nextPoint()) {
        nScanned++;
        //if(PDS.withCals)
        //        PDS.recalibrateEnergy();
        //if(PDS.fPID==PID_BETA && PDS.fType==TYPE_0_EVENT) {
        //        runCounts.add(PDS.getRun(),1.0);
        //        totalCounts[afp][gv]++;
        //}
        fillCoreHists(PDS,PDS.physicsWeight);
    }
    
    *runTimes += PDS.runTimes;
    *runCounts += PDS.runCounts;
    printf("\tscanned %i points\n",nScanned);
}

void RunAccumulator::makeOutput(bool doPlots) {
    printf("Generating output from %.0f counts over %.2f hours.\n", runCounts->GetTotal(), runTimes->GetTotal()/3600.);
    calculateResults();
    makeAnaResults();
    if(doPlots)
        makePlots();
    writeROOT();
}

unsigned int RunAccumulator::mergeDir(const string& d) {
    vector<string> inflnames;
    for(auto const& f: listdir(d)) inflnames.push_back(d+"/"+f+"/"+f+".root");
    size_t nmerged = addFiles(inflnames);
    makeOutput();
    return nmerged;
}

TH1* RunAccumulator::hToRate(TH1* h, int scaleAxes) {
    TH1* hc = (TH1*)h->Clone((h->GetName()+string("_Rate")).c_str());
    hc->Scale(1./runTimes->GetTotal());
    if(scaleAxes>=1) {
        for(int ix=1; ix<hc->GetNbinsX(); ix++) {
            double bx = hc->GetXaxis()->GetBinWidth(ix);
            if(scaleAxes>=2) {
                for(int iy=1; iy<hc->GetNbinsY(); iy++) {
                    Int_t b = hc->GetBin(ix,iy);
                    double by = hc->GetYaxis()->GetBinWidth(iy);
                    hc->SetBinContent(b, hc->GetBinContent(b)/bx/by);
                    hc->SetBinError(b, hc->GetBinError(b)/bx/by);
                }
            } else {
                hc->SetBinContent(ix, hc->GetBinContent(ix)/bx);
                hc->SetBinError(ix, hc->GetBinError(ix)/bx);
            }
        }
    }
    return hc;
}

/* --------------------------------------------------- */


void FGBGRegionsHist::setTemplate(const TH1& hTemplate, const RangeCutSet& rs) {
    myRegions = rs;
    for(size_t i=0; i<myRegions.regions.size(); i++)
        hs.push_back(myP->registerSavedHist(string(hTemplate.GetName())+"_"+to_str(i),hTemplate));
}

void FGBGRegionsHist::fill(double cutval, double x, double w) {
    for(size_t i=0; i<myRegions.regions.size(); i++) {
        for(auto const& cut: myRegions.regions[i])
            if(cut.first < cutval && cutval <= cut.second)
                hs[i]->Fill(x,w);
    }
}

void FGBGRegionsHist::fill(double cutval, double x, double y, double w) {
    for(size_t i=0; i<myRegions.regions.size(); i++) {
        for(auto const& cut: myRegions.regions[i])
            if(cut.first < cutval && cutval <= cut.second)
                dynamic_cast<TH2*>(hs[i])->Fill(x,y,w);
    }
}

void FGBGRegionsHist::makeRates(int axesScale, double xscale, bool bgsub) {
    hRates.resize(hs.size());
    for(size_t i=0; i<hs.size(); i++) {
        delete hRates[i];
        hRates[i] = myP->hToRate(hs[i],axesScale);
        hRates[i]->SetLineColor(hs.size()==2? 4-2*i : 1+i);
        hRates[i]->Scale(xscale/myRegions.norms[i]);
        if(bgsub && i) hRates[i]->Add(hRates[0],-1.);
    }
}

