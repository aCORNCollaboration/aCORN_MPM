#ifndef CONSTRUCT_SCINTILLATOR_HH
#define CONSTRUCT_SCINTILLATOR_HH

#include "Construction_Utils.hh"

/// nitrogen volume with scintillators
class ScintillatorConstruction: public MaterialUser {
public:
	/// constructor
	ScintillatorConstruction();
	
	G4double scint_Radius;			//< scintillator disc radius
	G4double scint_thick;			//< scintillator disc thickness

	G4LogicalVolume* scint_log;		//< scintillator logical volume
	
	/// construct logical container volume
	void Construct();
};

#endif
