#include "RunAccumulator.hh"
#include "GraphUtils.hh"
#include "SMExcept.hh"
#include "PathUtils.hh"
#include <time.h>

RunAccumulator::RunAccumulator(OutputManager* pnt, const std::string& nm, const std::string& inflName):
SegmentSaver(pnt,nm,inflName), isSimulated(false) {
    
    // initialize blind time to 0
    zeroCounters();
    
    // load existing data (if any)
    if(fIn) {
        QFile qOld(inflname+".txt");
        // fetch run counts, run times
        runCounts += TagCounter<RunID>(qOld.getFirst("runCounts"));
        runTimes += TagCounter<RunID>(qOld.getFirst("runTimes"));
    }
}

RunAccumulator::~RunAccumulator() {
    for(std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++)
        delete it->second;
}

void RunAccumulator::addPlugin(AnalyzerPlugin* AP) {
    if(myPlugins.find(AP->name) != myPlugins.end()) {
        SMExcept e("DuplicatePluginName");
        e.insert("myName",name);
        e.insert("pluginName",AP->name);
        throw(e);
    }
    myPlugins.insert(std::make_pair(AP->name,AP));
}

AnalyzerPlugin* RunAccumulator::getPlugin(const std::string& nm) {
    std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.find(nm);
    return it != myPlugins.end() ? it->second : NULL;
}

void RunAccumulator::fillCoreHists(BaseDataScanner& PDS, double weight) {
    for(std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++)
        it->second->fillCoreHists(PDS,weight);
}

void RunAccumulator::calculateResults() {
    printf("Calculating results for %s...\n",name.c_str());
    if(isCalculated) printf("*** Warning: repeat calculation!\n");
                            for(std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++) {
                                printf("... results in '%s' ...\n",it->first.c_str());
    it->second->calculateResults();
                            }
                            printf("Done calculating results for %s.\n",name.c_str());
        isCalculated = true;
}

void RunAccumulator::makePlots() {
    defaultCanvas->cd();
    if(!isCalculated) calculateResults();
                   printf("Generating plots for %s...\n",name.c_str());
    for(std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++) {
        printf("... plots in '%s' ...\n",it->first.c_str());
        it->second->makePlots();
    }
    printf("Done generating plots for %s...\n",name.c_str());
}

void RunAccumulator::compareMCtoData(RunAccumulator& OAdata) {
    defaultCanvas->cd();
    if(!isCalculated) calculateResults();
                                                    if(!OAdata.isCalculated) OAdata.calculateResults();
                                                    printf("Comparing MC %s and data %s...\n",name.c_str(),OAdata.name.c_str());
    for(std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++) {
        AnalyzerPlugin* AP = OAdata.getPlugin(it->second->name);
        if(AP) {
            printf("... comparison in '%s' ...\n",it->first.c_str());
            it->second->compareMCtoData(AP);
        }
    }
    printf("Done comparing MC %s and data %s.\n",name.c_str(),OAdata.name.c_str());
}

void RunAccumulator::zeroCounters() {
    runCounts = TagCounter<RunID>();
    runTimes = TagCounter<RunID>();
}

void RunAccumulator::addSegment(const SegmentSaver& S) {
    // histograms add
    SegmentSaver::addSegment(S);
    // recast
    const RunAccumulator& RA = dynamic_cast<const RunAccumulator&>(S);
    if(RA.isSimulated) isSimulated = true;
                                                    // add run counts, times
    runCounts += RA.runCounts;
    runTimes += RA.runTimes;
}

void RunAccumulator::scaleData(double s) {
    SegmentSaver::scaleData(s);
    runCounts.scale(s);
}

void RunAccumulator::write(std::string outName) {
    printf("Writing data to file '%s'...\n",outName.c_str());
    
    // clear previous tallies
    qOut.erase("runCounts");
    qOut.erase("runTimes");
    
    // record run counts, times
    qOut.insert("runCounts",runCounts.toStringmap());
    qOut.insert("runTimes",runTimes.toStringmap());
    
    // base class write
    SegmentSaver::write(outName);
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
    calculateResults();
    if(doPlots)
        makePlots();
    write();
    setWriteRoot(true);
}

unsigned int RunAccumulator::mergeDir() {
    std::vector<std::string> fnames = listdir(basePath);
    unsigned int nMerged = 0;
    for(std::vector<std::string>::iterator it = fnames.begin(); it != fnames.end(); it++) {
        // check whether data directory contains cloneable subdirectories
        std::string datinfl = basePath+"/"+(*it)+"/"+(*it);
        if(!inflExists(datinfl)) continue;
        SegmentSaver* subRA = makeAnalyzer(*it,datinfl);
        addSegment(*subRA);
        delete(subRA);
        nMerged++;
    }
    makeOutput();
    return nMerged;
}


TH1* RunAccumulator::hToRate(TH1* h, int scaleAxes) {
    TH1* hc = (TH1*)h->Clone();
    hc->Sumw2();
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

TH1* AnalyzerPlugin::registerHist(const std::string& hname, const std::string& title, unsigned int nbins, float xmin, float xmax) {
    myHists.push_back(myA->registerSavedHist(hname,title,nbins,xmin,xmax));
    return myHists.back();
}

TH1* AnalyzerPlugin::registerHist(const std::string& nm, const TH1& hTemplate) {
    myHists.push_back(myA->registerSavedHist(nm,hTemplate));
    return myHists.back();
}

/*------------------------------------------------------*/

FGBGRegionsHist::FGBGRegionsHist(AnalyzerPlugin* P): myP(P) {
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
    std::string hname = hTemplate.GetName();
    h[false] = myP->registerHist(hname+"_bg",hTemplate);
    h[true]  = myP->registerHist(hname+"_fg",hTemplate);
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

