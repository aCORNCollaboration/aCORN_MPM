/// \file PluginInterpolator.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "PluginInterpolator.hh"

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

