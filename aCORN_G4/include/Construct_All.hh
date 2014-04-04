//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.hh,v 1.11 2011-10-03 04:12:16 mmendenhall Exp $
// GEANT4 tag $Name:  $
//

#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "DetectorConstructionUtils.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "Rtypes.h"
#include <TString.h>
#include "TrackerSD.hh"
#include "Field.hh"

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
		
private:
	/// construct detector (Electro-)Magnetic Field
	void ConstructField(const TString filename);  
	
	// sensitive volumes
	TrackerSD* scint_SD[2];
	TrackerSD* Dscint_SD[2];
	TrackerSD* backing_SD[2];
	
	TrackerSD* winIn_SD[2];
	TrackerSD* winOut_SD[2];
	TrackerSD* trap_win_SD[2];
	TrackerSD* kevlar_SD[2];
	
	TrackerSD* mwpc_SD[2];
	TrackerSD* mwpc_planes_SD[2];
	TrackerSD* mwpcDead_SD[2];
	
	TrackerSD* source_SD;
	TrackerSD* trap_monitor_SD[2];
	
	TrackerSD* hall_SD;
	
	TrackerSD* siDet_SD;
	
	// UI commands
	G4UIdirectory* fDetectorDir;					//< UI Directory for detector-related commands
	
	G4UIcmdWithAString* fDetectorGeometry;			//< which detector geometry to construct
	G4String sGeometry;
	
	G4UIcmdWithAString* fFieldCmd;					//< whether to turn on/off the magnetic field
	G4String fieldSwitch;
	
	G4UIcmdWithABool* fAFPFieldCmd;					//< whether to enable the AFP fringe field
	bool fAddAFPField;
	
	G4UIcmdWithAString* fFieldMapFileCmd;			//< which field map to use
	TString sFieldMapFile;	
	
	G4UIcmdWith3VectorAndUnit* fSourceHolderPosCmd;	//< source holder position
	G4ThreeVector fSourceHolderPos;
	
	G4UIcmdWith3VectorAndUnit* fDetOffsetCmd;		//< Symmetrical detector offset from center axis
	G4ThreeVector fDetOffset;
	
	G4UIcmdWithADouble* fDetRotCmd;					//< Symmetrical detector rotation angle around Z axis (radians)
	Float_t fDetRot;
	
	G4UIcmdWithABool* fInFoilCmd;					//< construction of Indium 10um Al source foil
	bool makeInFoil;
	
	G4UIcmdWithADoubleAndUnit* fVacuumLevelCmd;		//< SCS bore vacuum
	Float_t fVacuumPressure;
	
	G4UIcmdWithADoubleAndUnit* fScintStepLimitCmd;	//< step size limiter in scintillator
	Float_t fScintStepLimit;
	
	G4UIcmdWithADoubleAndUnit* fMWPCBowingCmd;		//< additional width of MWPC due to window bowing
	Float_t fMWPCBowing;
	
	G4UIcmdWithADoubleAndUnit* fSourceFoilThickCmd;	//< source foil full thickness
	Float_t fSourceFoilThick;

	G4UIcmdWithADouble* fCrinkleAngleCmd;			//< decay trap foil crinkle angle
	Float_t fCrinkleAngle;
	
	G4UIcmdWithADouble* fMatterScaleCmd[2];			//< matter interaction scaling factor
	G4double fMatterScale[2];
	
	/// turn field on/off
	void SetFieldOnOff(G4String);
	
	Field* fpMagField;
};

#endif

