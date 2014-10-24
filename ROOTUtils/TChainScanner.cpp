#include "TChainScanner.hh"
#include "SMExcept.hh"
#include <stdlib.h>
#include <time.h>

TChainScanner::TChainScanner(const std::string& treeName): nEvents(0), nFiles(0), noEmpty(false), Tch(new TChain(treeName.c_str())),
currentEvent(0), noffset(0), nLocalEvents(0) {
    Tch->SetMaxVirtualSize(10000000);
}


int TChainScanner::addFile(const std::string& filename) {
    unsigned int oldEvents = nEvents;
    int nfAdded = Tch->Add(filename.c_str(),0);
    if(!nfAdded) {
        SMExcept e("missingFiles");
        e.insert("fileName",filename);
        throw e;
    }
    nEvents = Tch->GetEntries();
    nnEvents.push_back(nEvents-oldEvents);
    if(!nnEvents.back() && noEmpty) {
        SMExcept e("noEventsInFile");
        e.insert("fileName",filename);
        e.insert("nFiles",nfAdded);
        throw e;
    }
    if(!nFiles)
        setReadpoints(Tch);
    nFiles+=nfAdded;
    return nfAdded;
}

void TChainScanner::gotoEvent(unsigned int e) {
    currentEvent = e;
    Tch->GetEvent(currentEvent);
    smassert(Tch->GetTree());
    nLocalEvents = noffset = 0;
}

void TChainScanner::startScan(bool startRandom) { 
    if(!nEvents) {
        printf("Starting scan with no data... ");
        fflush(stdout);
        return;
    }
    if(startRandom) {
        if(!currentEvent) {
            srand(time(NULL));	// random random seed
            gotoEvent(rand()%Tch->GetEntries());
            printf("Scan Starting at offset %i/%i: ",currentEvent,nEvents);
        } else {
            printf("Scan Continuing at offset %i/%i: ",currentEvent,nEvents);
        }
    } else {
        gotoEvent(0);
        currentEvent = -1;
        printf(">%i< ",nEvents);
    }
    fflush(stdout);
}

void TChainScanner::SetBranchAddress(TTree* T, const std::string& bname, void* bdata) {
    smassert(bdata);
    smassert(T);
    Int_t err = T->SetBranchAddress(bname.c_str(),bdata);
    if(err && err != TTree::kNoCheck) {
        SMExcept e("TTreeBranchError");
        e.insert("branchName", bname);
        e.insert("errCode",err);
        throw e;
    }
}

void TChainScanner::speedload(unsigned int e, bool loadBaskets) {
    if(e < noffset || e-noffset >= nLocalEvents) {
        Tch->LoadTree(e);
        Tch->GetTree()->SetMaxVirtualSize(1000000);
        if(loadBaskets) Tch->GetTree()->LoadBaskets();
        nLocalEvents = Tch->GetTree()->GetEntries();
        noffset = Tch->GetChainOffset();
    }
    Tch->GetTree()->GetEvent(e-noffset);
}

bool TChainScanner::nextPoint() {
    if(!nEvents) return false;
    ++currentEvent;
    if(currentEvent >= nEvents) {
        printf("\n");
        startScan();
        return false;
    }
    if(!(currentEvent%(nEvents/20))) {
        printf("*"); fflush(stdout);
    }
    speedload(currentEvent);
    return true;
}
