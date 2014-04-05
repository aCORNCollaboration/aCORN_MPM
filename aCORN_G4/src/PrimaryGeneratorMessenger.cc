#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* myGun): Action(myGun) {
	gunDir = new G4UIdirectory("/benchmark/gun/");
	gunDir->SetGuidance("PrimaryGenerator control");
	
	gunTypeCmd = new G4UIcmdWithAString("/benchmark/gun/type",this);
	gunTypeCmd->SetGuidance("Set the generator gun type.");
	gunTypeCmd->SetGuidance(" Choices: eGun eGunRandMomentum");
	gunTypeCmd->SetDefaultValue("eGun");
	gunTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	gunPtclCmd = new G4UIcmdWithAString("/benchmark/gun/particle",this);
	gunPtclCmd->SetGuidance("Set the gun particle thrown.");
	gunPtclCmd->SetDefaultValue("e-");
	gunPtclCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	positionerCmd = new G4UIcmdWithAString("/benchmark/gun/positioner",this);
	positionerCmd->SetGuidance("Set the generator gun positioner.");
	positionerCmd->SetGuidance(" Choice : Fixed");
	positionerCmd->SetDefaultValue("Fixed");
	positionerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	eventFileCmd = new G4UIcmdWithAString("/benchmark/gun/evtfile",this);
	eventFileCmd->SetGuidance("Set input file for events");
	eventFileCmd->SetDefaultValue("");
	eventFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}


PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger() {
	delete gunTypeCmd;
	delete positionerCmd;
	delete gunDir;
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 
	if( command == gunTypeCmd )
		Action->SetGunType(newValue);
	if( command == gunPtclCmd )
		Action->SetParticleType(newValue);
	if( command == positionerCmd )
		Action->SetPositioner(newValue);
	if( command == eventFileCmd )
		Action->SetEventFile(newValue);
}
