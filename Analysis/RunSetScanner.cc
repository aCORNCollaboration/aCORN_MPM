#include "RunSetScanner.hh"
#include "PathUtils.hh"
#include "SMExcept.hh"
#include <stdio.h>
#include <stdlib.h>

RunSetScanner::RunSetScanner(const std::string& treeName):
TChainScanner(treeName), totalTime(0) { }

RunSetScanner::~RunSetScanner() { }

unsigned int RunSetScanner::addRuns(const std::vector<RunNum>& rns) {
    printf("\n------------------- Assembling %i runs into TChain... ",(int)rns.size()); fflush(stdout);
    unsigned int n = 0;
    for(std::vector<RunNum>::const_iterator it = rns.begin(); it != rns.end(); it++) {
        n+=addRun(*it);
        printf("*"); fflush(stdout);
    }
    printf("------------------- %i Runs, %i events found, %.2fh running.\n",getnFiles(),nEvents,totalTime/3600.0);
    return n;
}

void RunSetScanner::display() {
    printf("RunSetScanner: %i runs, %i events\n",getnFiles(),nEvents);
}

void RunSetScanner::speedload(unsigned int e) {
    if(e < noffset || e-noffset >= nLocalEvents) {
        Tch->LoadTree(e);
        Tch->GetTree()->LoadBaskets();
        nLocalEvents = Tch->GetTree()->GetEntries();
        noffset = Tch->GetChainOffset();
        if((int)runlist.size()>Tch->GetTreeNumber())
            evtRun = runlist[Tch->GetTreeNumber()];
        else
            evtRun = Tch->GetTreeNumber();
        loadNewRun(evtRun);
    }
    Tch->GetTree()->GetEvent(e-noffset);
}

bool RunSetScanner::addRun(RunNum r) {
    std::string f = locateRun(r);
    if(f.size() && addFile(f)) {
        runlist.push_back(r);
        //RunInfo R = CalDBSQL::getCDB()->getRunInfo(r);
        //BlindTime b = CalDBSQL::getCDB()->fiducialTime(r);
        //totalTime += b;
        runTimes.add(r,0);
        return true;
    }
    printf("**** FAILED TO LOCATE analyzed data for run %i at '%s'! *****\n",r,f.c_str());
    return false;
}
