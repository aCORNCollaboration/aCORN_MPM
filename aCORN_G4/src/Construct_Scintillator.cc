#include "Construct_Scintillator.hh"
#include "G4SystemOfUnits.hh"
#include <G4Polycone.hh>
#include <cassert>

ScintillatorConstruction::ScintillatorConstruction(): scint_Radius(7.5*cm), scint_thick(8.0*mm) { }

void ScintillatorConstruction::Construct() {
	
	G4Tubs* scint_tube = new G4Tubs("scint_tube",	// name
									0.,				// rmin
									scint_Radius,	// rmax
									scint_thick/2.,	// half z
									0.,				// start angle
									2*M_PI			// total angle
									);
	scint_log = new G4LogicalVolume(scint_tube, Sci, "scint_log");
	G4VisAttributes* visScint= new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.2));
	scint_log->SetVisAttributes(visScint);
}
