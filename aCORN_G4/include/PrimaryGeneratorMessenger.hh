#ifndef PRIMARYGENERATORMESSENGER_HH
#define PRIMARYGENERATORMESSENGER_HH

#include "PrimaryGeneratorAction.hh"

#include <globals.hh>
#include <G4UImessenger.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIdirectory.hh>

class PrimaryGeneratorAction;

/// messenger UI for primary event generator
class PrimaryGeneratorMessenger: public G4UImessenger {
public:
	/// constructor
    PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
	/// destructor
	~PrimaryGeneratorMessenger();
    
	/// receive command
    void SetNewValue(G4UIcommand*, G4String);
    
private:
    PrimaryGeneratorAction* Action;				//< primary generator action being controlled
    G4UIdirectory*			gunDir;				//< '/gun/' commands directory
    G4UIcmdWithAString*		gunTypeCmd;			//< control momentum generator (beam, random, file-specified)
	G4UIcmdWithAString*		gunPtclCmd;			//< control default particle to throw
	G4UIcmdWithAString*		positionerCmd;		//< control event positioner to use
	G4UIcmdWithAString*		eventFileCmd;		//< control file to read events from
};

#endif
