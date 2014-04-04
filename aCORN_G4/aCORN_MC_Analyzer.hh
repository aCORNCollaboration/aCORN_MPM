#include "AnalyzerBase.hh"
#include "Enums.hh"
#include "WirechamberCalibrator.hh"

/// number of sensitive detector regions
#define N_SD 1

class aCORN_MC_Analyzer: public aCORN_G4_Analyzer {
public:
	/// constructor
	aCORN_MC_Analyzer(const std::string& outfname);
	
	bool saveAllEvents;						//< whether to save non-energy-depositing events to file
	
	Double_t Edep;							//< scintillator deposited energy
	Double_t EdepQ;							//< quenched energy in scintillator
	Double_t ScintPos[Z_DIRECTION+1];		//< scintillator deposited energy weighted position
	Double_t ScintPosSigma[Z_DIRECTION+1];	//< scintillator quenched energy weighted position variance
	
	Double_t EdepSD[N_SD];					//< array for energy deposition in all SDs
	Double_t thetaInSD[N_SD];				//< entrance angle in each sensitive detector
	Double_t thetaOutSD[N_SD];				//< exit angle for each sensitive detector
	Double_t keInSD[N_SD];					//< kinetic energy entering each sensitive detector
	Double_t keOutSD[N_SD];					//< kinetic energy exiting each sensitive detector
	Int_t hitCountSD[N_SD];					//< count of primary tracks in each volume
	Double_t EdepAll;						//< total edep in all SDs
	Double_t hitTimeSD[N_SD];				//< earliest hit time in each SD
		
protected:
	
	/// add additional branches to output tree
	virtual void setupOutputTree();
	
	/// reset analysis values for new event
	virtual void resetAnaEvt();
	/// process current track segment
	virtual void processTrack();
	/// final whole-event processing
	virtual void processEvent();
	/// determine whether an event should be saved to output file
	virtual bool saveEvent() { return saveAllEvents || Edep > 0; }
};
