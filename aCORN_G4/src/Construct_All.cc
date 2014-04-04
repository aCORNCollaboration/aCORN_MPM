#include "Construct_All.hh"
#include "G4SystemOfUnits.hh"

#include "AnalysisManager.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"
#include "G4TransportationManager.hh"

#include "G4UserLimits.hh"
#include "G4PVParameterised.hh"

#include "TrackerSD.hh"
#include "Field.hh"

#include <cassert>

DetectorConstruction::DetectorConstruction(): experimentalHall_log(NULL), experimentalHall_phys(NULL), fpMagField(NULL) {
	
	fDetectorDir = new G4UIdirectory("/detector/");
	fDetectorDir->SetGuidance("/detector control");
	
	fFieldMapFileCmd = new G4UIcmdWithAString("/detector/fieldmapfile",this);
	fFieldMapFileCmd->SetGuidance("Set B field map file");
	sFieldMapFile = "";
	
	fVacuumLevelCmd = new G4UIcmdWithADoubleAndUnit("/detector/vacuum",this);
	fVacuumLevelCmd->SetGuidance("Set SCS vacuum pressure");
	fVacuumPressure = 0;
	
	fScintStepLimitCmd = new G4UIcmdWithADoubleAndUnit("/detector/scintstepsize",this);
	fScintStepLimitCmd->SetGuidance("step size limit in scintillator");
	fScintStepLimitCmd->SetDefaultValue(0.1*mm);

}

void DetectorConstruction::SetNewValue(G4UIcommand * command, G4String newValue) {
	if (command == fFieldMapFileCmd) {
		sFieldMapFile = newValue;
	} else if (command == fVacuumLevelCmd) {
		fVacuumPressure = fVacuumLevelCmd->GetNewDoubleValue(newValue);
	} else if (command == fScintStepLimitCmd) {
		fScintStepLimit = fScintStepLimitCmd->GetNewDoubleValue(newValue);
		G4cout << "Setting step limit in solids to " << fScintStepLimit/mm << "mm" << G4endl;
	} else {
		G4cerr << "Unknown command:" << command->GetCommandName() << " passed to DetectorConstruction::SetNewValue\n";
    }
}

TrackerSD* registerSD(G4String sdName) {
	TrackerSD* sd = new TrackerSD(sdName);
	G4SDManager::GetSDMpointer()->AddNewDetector(sd);
	gAnalysisManager->SaveSDName(sdName);
	return sd;
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
	
	//////////////////////////////////////
	// init materials
	//////////////////////////////////////
	setVacuumPressure(fVacuumPressure);
	
	//////////////////////////////////////
	// user step limits
	//////////////////////////////////////
	G4UserLimits* myCoarseLimits = new G4UserLimits();
	myCoarseLimits->SetMaxAllowedStep(10*m);
	G4UserLimits* mySolidLimits = new G4UserLimits();
	mySolidLimits->SetMaxAllowedStep(fScintStepLimit);
		
	///////////////////////////////////////
	//experimental Hall
	///////////////////////////////////////
	const G4double expHall_halfx=1.0*m;
	const G4double expHall_halfy=1.0*m;
	const G4double expHall_halfz=4.0*m;  	
	G4Box* experimentalHall_box = new G4Box("expHall_box",expHall_halfx,expHall_halfy,expHall_halfz);
	experimentalHall_log = new G4LogicalVolume(experimentalHall_box,Vacuum,"World_Log");  
	experimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);
	experimentalHall_log->SetUserLimits(myCoarseLimits);
	experimentalHall_phys = new G4PVPlacement(NULL, G4ThreeVector(), "World_Phys", experimentalHall_log, NULL, false, 0);
	
	
	////////////////////////////////////////
	// detector components
	////////////////////////////////////////
	scint.Construct();
	scint_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0,-0.3*m), "scint_phys", scint.scint_log, experimentalHall_phys, false, 0);
		
	////////////////////////////////////////
	// sensitive volumes
	////////////////////////////////////////
	scint_SD = registerSD("scint_SD");
	scint.scint_log->SetSensitiveDetector(scint_SD);
			
	
	// construct magnetic field
	cout<<"##### "<<sFieldMapFile<<" #####"<<endl;
	ConstructField(sFieldMapFile);
	
	return experimentalHall_phys;
}

#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4SimpleHeum.hh"
#include "G4HelixHeum.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4HelixMixedStepper.hh"

void DetectorConstruction::ConstructField(const G4String& filename) {
	
	static G4bool fieldIsInitialized = false;
	
	if(!fieldIsInitialized) {
		cout<<"##### Constructing Field #####"<<endl;
		
		// get magnetic field profile
		fpMagField = new Field(filename);
		// set up field manager for this profile
		G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
		fieldMgr->SetDetectorField(fpMagField);
		// set up default chord finder
		fieldMgr->CreateChordFinder(fpMagField);
		
		// Select stepper
		G4MagIntegratorStepper* pStepper;
		G4Mag_UsualEqRhs* fEquation = new G4Mag_UsualEqRhs(fpMagField); // equation of motion in magnetic field
		//pStepper = new G4ClassicalRK4 (fEquation);		// general case for "smooth" EM fields
		//pStepper = new G4SimpleHeum( fEquation );			// for slightly less smooth EM fields
		//pStepper = new G4HelixHeum( fEquation );			// for "smooth" pure-B fields
		//pStepper = new G4HelixImplicitEuler( fEquation );	// for less smooth pure-B fields; appears ~50% faster than above
		//pStepper = new G4HelixSimpleRunge( fEquation );	// similar speed to above
		//pStepper = new G4HelixExplicitEuler( fEquation );	// about twice as fast as above
		pStepper = new G4HelixMixedStepper(fEquation,6);	// avoids "Stepsize underflow in Stepper" errors
		fieldMgr->GetChordFinder()->GetIntegrationDriver()->RenewStepperAndAdjust(pStepper);
		
		// set required accuracy for finding intersections
		fieldMgr->GetChordFinder()->SetDeltaChord(100.0*um);
		// set integration relative error limits for small and large steps
		fieldMgr->SetMinimumEpsilonStep(1e-6);
		fieldMgr->SetMaximumEpsilonStep(1e-5);
		// set integration absolute error limit
		fieldMgr->SetDeltaOneStep(0.1*um);
		// allow lots of looping
		G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetMaxLoopCount(INT_MAX);
		
		fieldIsInitialized = true;
	}
}
