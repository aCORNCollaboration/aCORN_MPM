#ifndef ANALYZERBASE_HH
#define ANALYZERBASE_HH

#include <vector>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cassert>

#include "ROOT_Headers.hh"

#include "MCEvent.hh"
#include "TrackInfo.hh"
#include "PrimaryInfo.hh"

using namespace std;

/// base class for analyzing aCORN Geant4 MC simulations
class aCORN_G4_Analyzer {
public:
	/// constructor
	aCORN_G4_Analyzer(const std::string& outfname);
	/// destructor
	virtual ~aCORN_G4_Analyzer();
	
	/// analyze a file
	void analyzeFile(const string& fname);
	/// analyze all files in list file
	void analyzeFileList(const string& flist);
	
	Int_t fTrapped;				//< whether this particle is trapped in mag field (lives too long)
	Double_t fCompTime;			//< computer time needed to track this event
	Long_t seed;				//< random seed used for this event
	
	Double_t primKE;			//< primary event kinetic energy
	Double_t primTheta;			//< primary event emission angle
	Double_t primWeight;		//< primary event generator weight
	Double_t primPos[4];		//< primary event position (4th coordinate = radius)
	int pID;					//< track pID
	int detectorID;				//< track detector ID
	int trackID;				//< track segment number
	
protected:
	
	TTree* anaTree;				//< analysis results output tree
	TFile* outf;				//< output file
	MCEvent* myevt;			//< current event being analyzed
	TrackInfo* trackinfo;		//< current track info
	PrimaryInfo* priminfo;	//< current primary info
	
	/// add additional branches to output tree
	virtual void setupOutputTree() {}
	
	/// reset analysis values for new event
	virtual void resetAnaEvt() {}
	/// process current track segment
	virtual void processTrack() {}
	/// final whole-event processing
	virtual void processEvent() {}
	/// determine whether an event should be saved to output file
	virtual bool saveEvent() { return true; }
};

#endif
