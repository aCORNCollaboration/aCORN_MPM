#ifndef CONSTRUCT_SCINTILLATOR_HH
#define CONSTRUCT_SCINTILLATOR_HH 1

#include "Construction_Utils.hh"

/// nitrogen volume with scintillators
class ScintillatorConstruction: public MaterialUser {
public:
	/// constructor
	ScintillatorConstruction(): scint_Radius(7.5*cm), scint_thick(8.0*mm) { }
	
	G4double scint_Radius;				//< scintillator disc radius
	G4double scint_thick;				//< scintillator disc thickness

	G4LogicalVolume* scint_log;			//< scintillator logical volume
	
	/// construct logical container volume
	void Construct();
	
protected:
	G4VPhysicalVolume* scint_phys;		//< scintillator physical volume
};

#endif
