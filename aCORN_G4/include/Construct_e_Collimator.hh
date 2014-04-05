#ifndef CONSTRUCT_E_COLLIMATOR_HH
#define CONSTRUCT_E_COLLIMATOR_HH

#include "Construction_Utils.hh"
#include <vector>

/// nitrogen volume with scintillators
class eCollimatorConstruction: public MaterialUser {
public:
	/// constructor
	eCollimatorConstruction();
	
	G4double r_inner;			//< inner radius
	G4double r_outer;			//< outer radius
	
	G4LogicalVolume* eCollimator_log;		//< overall logical volume
	
	/// construct logical container volume
	void Construct();
	/// get constructed unit length
	G4double getLength() const { return length; }
	
protected:
	
	G4double length;			//< length
	G4double c_thick;			//< ring thickness
	
	G4LogicalVolume* aperture_log;				//< logical volume for one collimator ring
	std::vector<G4VPhysicalVolume*> apertures;	//< tungsten collimator rings
};

#endif
