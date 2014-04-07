#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH

#include "ROOT_Headers.hh"

class G4VPhysicalVolume;
class G4Event;
class G4Run;
class G4Track;
class G4Step;
class G4PrimaryVertex;
class G4PrimaryParticle;

#include <globals.hh>

#include "EventAction.hh"
#include "SteppingAction.hh"
#include "TrackerHit.hh"
#include "TrackerSD.hh"
#include "MCEvent.hh"
#include <vector>

class AnalysisManager;

extern AnalysisManager *gAnalysisManager; // global AnalysisManager

/// handles storing ROOT data for each event
class AnalysisManager {
	
public:
	/// constructor
	AnalysisManager();
	/// destructor
	~AnalysisManager() {
		if (gAnalysisManager == this)
			gAnalysisManager = (AnalysisManager *)0;
	}
	/// get global analysis manager
	static AnalysisManager* GetAnalysisManager() { return gAnalysisManager; }
	
	/// open ROOT output file
	void OpenFile(const G4String filename);
	/// write and close ROOT output file
	void CloseFile();
	/// clear event info
	void Clear(){ mcEvent.ClearEvent(); }
	/// set up output tree
	void CreateTrees();
	
	void StoreHitCollectionIDs();
	
	void SetRunNumber(const Int_t runno) { fRunNumber = runno; }
	Int_t GetRunNumber(void) const {return fRunNumber;}
	
	/// convert primary events data to ROOT form
	void FillPrimaryData(const G4Event* evt_in, const long);
	/// convert tracking data to ROOT form
	void FillTrackerData(const G4Event *evt);
	/// fill output tree
	void FillEventTree();
	
	MCEvent mcEvent;		//< saveable event
	MCEvent* pMcEvent;	//< pointer to said event for TTree branch setup
	
	/// store a sensitive detector name to record tracking info for
	void SaveSDName(const G4String name){ fSDNames.push_back(name); }
		
private:
	
	TFile *fROOTOutputFile;		//< ROOT output file
	TTree *fEventTree;			//< ROOT output TTree
	Int_t fRunNumber;  			//< MC run number
	vector<G4int> detectorIDs;	//< list of SD ID numbers
	vector<G4String> fSDNames;	//< list of SD names corresponding to ID numbers
};

#endif
