/// \file RunAccumulator.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "RunAccumulator.hh"
#include "GraphUtils.hh"
#include "SMExcept.hh"
#include "PathUtils.hh"
#include <time.h>
#include <TH1.h>

RunAccumulator::RunAccumulator(OutputManager* pnt, const std::string& nm, const std::string& inflName):
PluginSaver(pnt,nm,inflName), isSimulated(false) {
    TH1::SetDefaultSumw2(true); // all histograms default to having errorbars
    
    // initialize blind time to 0
    zeroCounters();
    
    // load existing data (if any)
    if(fIn) {
        SMFile qOld(dropLast(inflname,".")+".txt");
        // fetch run counts, run times
        runCounts += TagCounter<RunID>(qOld.getFirst("runCounts"));
        runTimes += TagCounter<RunID>(qOld.getFirst("runTimes"));
    }
}

void RunAccumulator::fillCoreHists(BaseDataScanner& PDS, double weight) {
    for(auto it = myBuilders.begin(); it != myBuilders.end(); it++) {
        auto P = dynamic_cast<RunAccumulatorPlugin*>(it->second->thePlugin);
        if(P) P->fillCoreHists(PDS,weight);
    }
}

AnaResult RunAccumulator::makeBaseResult() const {
    AnaResult baseResult;
    if(runCounts.counts.size()) {
        baseResult.start = runCounts.counts.begin()->first;
        baseResult.end = runCounts.counts.rbegin()->first;
    }
    return baseResult;
}

void RunAccumulator::makeAnaResults() {
    if(!isCalculated) calculateResults();
    AcornDB::ADB().beginTransaction();
    
    AnaResult baseResult = makeBaseResult();
    baseResult.value = runCounts.total();
    AcornDB::ADB().uploadAnaResult("total_counts", "Total analyzed counts", baseResult);
    baseResult.value = runTimes.total();
    AcornDB::ADB().uploadAnaResult("total_time", "Total analyzed time [s]", baseResult);
    //for(map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++) {
    //    it->second->makeAnaResults();
    //}
    assert(false);
    AcornDB::ADB().endTransaction();
}

void RunAccumulator::zeroCounters() {
    runCounts = TagCounter<RunID>();
    runTimes = TagCounter<RunID>();
}

void RunAccumulator::addSegment(const SegmentSaver& S) {
    // histograms add
    PluginSaver::addSegment(S);
    
    // recast
    const RunAccumulator& RA = dynamic_cast<const RunAccumulator&>(S);
    if(RA.isSimulated) isSimulated = true;
    // add run counts, times
    runCounts += RA.runCounts;
    runTimes += RA.runTimes;
}

void RunAccumulator::scaleData(double s) {
    PluginSaver::scaleData(s);
    runCounts.scale(s);
}

void RunAccumulator::write(string outName) {
    printf("Writing data to file '%s'...\n",outName.c_str());
    
    // clear previous tallies
    qOut.erase("runCounts");
    qOut.erase("runTimes");
    
    // record run counts, times
    qOut.insert("runCounts",runCounts.toStringmap());
    qOut.insert("runTimes",runTimes.toStringmap());
    
    if(!outName.size()) outName = basePath+"/"+name+".txt";
    qOut.commit(outName);
}

void RunAccumulator::loadProcessedData(BaseDataScanner& PDS) {
    printf("Loading processed data...\n");
    if(!PDS.getnFiles())
        return;
    PDS.startScan();
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
    
    runTimes += PDS.runTimes;
    runCounts += PDS.runCounts;
    printf("\tscanned %i points\n",nScanned);
}

void RunAccumulator::makeOutput(bool doPlots) {
    printf("Generating output from %.0f counts over %.2f hours.\n", runCounts.total(), runTimes.total()/3600.);
    calculateResults();
    makeAnaResults();
    if(doPlots)
        makePlots();
    write();
    writeROOT();
}

unsigned int RunAccumulator::mergeDir() {
    vector<string> fnames = listdir(basePath);
    unsigned int nMerged = 0;
    for(vector<string>::iterator it = fnames.begin(); it != fnames.end(); it++) {
        // check whether data directory contains cloneable subdirectories
        string datinfl = basePath+"/"+(*it)+"/"+(*it);
        if(!fileExists(datinfl)) continue;
        SegmentSaver* subRA = makeAnalyzer(*it,datinfl);
        addSegment(*subRA);
        delete(subRA);
        nMerged++;
    }
    makeOutput();
    return nMerged;
}

TH1* RunAccumulator::hToRate(TH1* h, int scaleAxes) {
    TH1* hc = (TH1*)h->Clone((h->GetName()+string("_Rate")).c_str());
    hc->Scale(1./runTimes.total());
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

FGBGRegionsHist::FGBGRegionsHist(RunAccumulatorPlugin* P): myP(P) {
    for(int i=0; i<2; i++) {
        h[i] = hRates[i] = NULL;
        totalLength[i] = 0;
    }
}

FGBGRegionsHist::~FGBGRegionsHist() {
    for(int i=0; i<2; i++) {
        if(hRates[i]) delete hRates[i];
        hRates[i] = NULL;
    }
}

void FGBGRegionsHist::setTemplate(const TH1& hTemplate) {
    string hname = hTemplate.GetName();
    h[false] = myP->registerSavedHist(hname+"_bg",hTemplate);
    h[true]  = myP->registerSavedHist(hname+"_fg",hTemplate);
}

void FGBGRegionsHist::addRegion(double x0, double x1, bool fg) {
    regions[fg].push_back(pair<double,double>(x0,x1));
    totalLength[fg] += x1-x0;
}

void FGBGRegionsHist::fill(double cutval, double x, double w) {
    for(int i=0; i<2; i++) {
        if(!h[i]) continue;
        for(auto it = regions[i].begin(); it != regions[i].end(); it++)
            if(it->first < cutval && cutval <= it->second)
                h[i]->Fill(x,w);
    }
}

void FGBGRegionsHist::fill(double cutval, double x, double y, double w) {
    for(int i=0; i<2; i++) {
        if(!h[i]) continue;
        for(auto it = regions[i].begin(); it != regions[i].end(); it++)
            if(it->first < cutval && cutval <= it->second)
                dynamic_cast<TH2*>(h[i])->Fill(x,y,w);
    }
}

void FGBGRegionsHist::makeRates(int axesScale, double xscale) {
    for(int i=0; i<2; i++) {
        if(hRates[i]) delete hRates[i];
        hRates[i] = myP->myA->hToRate(h[i],axesScale);
        hRates[i]->SetLineColor(4-2*i);
    }
    hRates[false]->Scale(xscale*totalLength[true]/totalLength[false]);
    hRates[true]->Scale(xscale);
    hRates[true]->Add(hRates[false],-1.);
}

