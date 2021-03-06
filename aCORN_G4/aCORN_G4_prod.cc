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
// $Id: ucnG4_prod.cc,v 1.9 2011-10-12 19:31:16 mmendenhall Exp $
// GEANT4 tag $Name:  $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Construct_All.hh"
#include "PhysicsList_495.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "SteppingAction_Verbose.hh"
#include "AnalysisManager.hh"

#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>
#include <G4RunManager.hh>
#include <G4UImanager.hh>
#include <G4UIExecutive.hh>
#ifdef G4VIS_USE
#include <G4VisExecutive.hh>
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv) {
	
	if(argc < 2) {
		G4cout << "Usage:" << G4endl << "\t" << argv[0] << " <macro filename>" << G4endl;
		return 0;
	}
	
	// User Verbose stepping output class
	G4VSteppingVerbose::SetInstance(new SteppingAction_Verbose());
	
	// Run manager
	G4RunManager* runManager = new G4RunManager;
	
	// User Initialization classes (mandatory)
	DetectorConstruction* detector = new DetectorConstruction();
	runManager->SetUserInitialization(detector);
	runManager->SetUserInitialization(new PhysicsList_495(false));
	
	new G4UnitDefinition("torr","torr","Pressure",atmosphere/760.);
	
#ifdef G4VIS_USE
	// Visualization, if you choose to have it!
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();
#endif
	
	// User Action classes
	runManager->SetUserAction(new PrimaryGeneratorAction(detector));
	runManager->SetUserAction(new RunAction);
	runManager->SetUserAction(new EventAction);
	runManager->SetUserAction(new SteppingAction);
	
	//create global analysis manager for histograms and trees
	gAnalysisManager = new AnalysisManager();
	
	// Execute input macro file
	G4UImanager * UI = G4UImanager::GetUIpointer(); 
	G4String command = "/control/execute ";
	G4String fileName = argv[1];
	UI->ApplyCommand(command+fileName);
	
	// interactive UI session
	if(argc >= 3 && std::string(argv[argc-1]) == "ui") {
		G4UIExecutive* UIuser = new G4UIExecutive(argc, argv);
		UIuser->SessionStart();
		delete UIuser;
	}
	
#ifdef G4VIS_USE
	delete visManager;
#endif
	delete runManager;
	
	return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

