#include "RunSetScanner.hh"
#include "PathUtils.hh"

#include <iostream>
#include <stdlib.h>

RunSetScanner::RunSetScanner(const std::string& treeName):
TChainScanner(treeName), totalTime(0) { }

RunSetScanner::~RunSetScanner() { }

unsigned int RunSetScanner::addRuns(const vector<RunID>& rns) {
    printf("\n------------------- Assembling %i runs into TChain... ",(int)rns.size()); fflush(stdout);
    unsigned int n = 0;
    for(vector<RunID>::const_iterator it = rns.begin(); it != rns.end(); it++) {
        n+=addRun(*it);
        printf("*"); fflush(stdout);
    }
    printf("------------------- %i Runs, %i events found, %.2fh running.\n",getnFiles(),nEvents,totalTime/3600.0);
    return n;
}

void RunSetScanner::display() {
    printf("RunSetScanner: %i runs, %i events\n",getnFiles(),nEvents);
}

void RunSetScanner::nextTreeLoaded() {
    assert((size_t)Tch->GetTreeNumber() < runlist.size());
    evtRun = runlist[Tch->GetTreeNumber()];
    loadNewRun(evtRun);
}

void RunSetScanner::speedload(unsigned int e) {
    TChainScanner::speedload(e);
    calibrate();
}

bool RunSetScanner::addRun(RunID r) {
    Tch->SetMaxVirtualSize(100000); // need this to prevent bad_alloc error on big files
    string f = locateRun(r);
    if(f.size() && addFile(f)) {
        runlist.push_back(r);
        double rtime = _getRunTime(r);
        totalTime += rtime;
        runTimes.add(r,rtime);
        runCounts.add(r,nnEvents.back());
        return true;
    }
    std::cout << "**** FAILED TO LOCATE analyzed data for run " << r << " at '" << f << "'! *****\n";
    return false;
}

