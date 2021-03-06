#ifndef CONSTRUCT_ALL_HH
#define CONSTRUCT_ALL_HH

#include "Construction_Utils.hh"
#include "Construct_Scintillator.hh"
#include "Construct_e_Collimator.hh"
#include "TrackerSD.hh"
#include "Field.hh"

#include <G4VUserDetectorConstruction.hh>
#include <G4UImessenger.hh>
#include <G4UIdirectory.hh>
#include <G4UIcommand.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UIcmdWith3VectorAndUnit.hh>
#include <G4UIcmdWithAString.hh>

class DetectorConstruction : public G4VUserDetectorConstruction, G4UImessenger, MaterialUser {
public:
	/// constructor
	DetectorConstruction();
	
	/// construct detector geometry
	G4VPhysicalVolume* Construct();
	/// UI interface
	virtual void SetNewValue(G4UIcommand * command,G4String newValue);
	
	// world volume
	G4LogicalVolume* experimentalHall_log;	
	G4VPhysicalVolume* experimentalHall_phys;
	
	// components
	G4VPhysicalVolume* scint_phys;		//< electron scintillator
	G4VPhysicalVolume* eCollim_phys;	//< electron collimator
	
private:

	ScintillatorConstruction scint;		//< electron scintillator
	eCollimatorConstruction eCollim;	//< electron collimator
	
	/// construct detector (Electro-)Magnetic Field
	void ConstructField(const G4String& filename);
	Field* fpMagField;
		
	// sensitive volumes
	TrackerSD* scint_SD;
	
	// UI commands
	G4UIdirectory* fDetectorDir;					//< UI Directory for detector-related commands
	
	G4UIcmdWithAString* fFieldMapFileCmd;			//< which field map to use
	G4String sFieldMapFile;
	
	G4UIcmdWithADoubleAndUnit* fVacuumLevelCmd;		//< apparatus vacuum pressure
	Float_t fVacuumPressure;
	
	G4UIcmdWithADoubleAndUnit* fScintStepLimitCmd;	//< step size limiter in scintillator
	Float_t fScintStepLimit;
};

#endif
