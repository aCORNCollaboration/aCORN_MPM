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

void RunSetScanner::speedload(unsigned int e, bool loadBaskets) {
    if(e < noffset || e-noffset >= nLocalEvents) {
        if(loadBaskets) Tch->GetTree()->DropBaskets();
        Tch->GetEvent(e);
        if(loadBaskets) Tch->GetTree()->LoadBaskets();
        nLocalEvents = Tch->GetTree()->GetEntries();
        noffset = Tch->GetChainOffset();
        assert((size_t)Tch->GetTreeNumber() < runlist.size());
        evtRun = runlist[Tch->GetTreeNumber()];
        loadNewRun(evtRun);
    } else {
        Tch->GetEvent(e);
    }
    // Tch->GetTree()->GetEvent(e-noffset);
    calibrate();
}

bool RunSetScanner::addRun(RunID r) {
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

