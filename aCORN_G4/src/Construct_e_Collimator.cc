#include "Construct_e_Collimator.hh"
#include "G4SystemOfUnits.hh"
#include <G4Polycone.hh>
#include <cassert>
#include "strutils.hh"

///////
// see Noid thesis, p. 78...
///////

eCollimatorConstruction::eCollimatorConstruction(): r_inner(5.5*cm), r_outer(10.0*cm), c_thick(0.05*cm) { }

void eCollimatorConstruction::Construct() {
	
	// spacing between 17 subsequent collimator aperture top surfaces
	// starting from the top, moving down to "post hole digger"
	std::vector<G4double> aperture_spacing;
	aperture_spacing.push_back(0*cm);
	// spacer ring
	aperture_spacing.push_back(1.025*cm);
	aperture_spacing.push_back(1.281*cm);
	aperture_spacing.push_back(1.538*cm);
	aperture_spacing.push_back(1.794*cm);
	aperture_spacing.push_back(2.050*cm);
	aperture_spacing.push_back(2.306*cm);
	aperture_spacing.push_back(2.563*cm);
	// spacer ring
	aperture_spacing.push_back(2.819*cm);
	aperture_spacing.push_back(3.075*cm);
	aperture_spacing.push_back(3.331*cm);
	aperture_spacing.push_back(3.588*cm);
	aperture_spacing.push_back(3.844*cm);
	// spacer ring
	aperture_spacing.push_back(4.100*cm);
	aperture_spacing.push_back(4.356*cm);
	aperture_spacing.push_back(4.613*cm);
	// spacer ring
	aperture_spacing.push_back(4.869*cm);
	
	// convert to cumulative position; get total length
	for(unsigned int i=1; i<aperture_spacing.size(); i++) aperture_spacing[i] += aperture_spacing[i-1];
	length = aperture_spacing.back() + c_thick;
	
	// overall volume
	G4Tubs* e_collim_tube = new G4Tubs("e_collim_tube",	// name
									   r_inner,			// rmin
									   r_outer,			// rmax
									   length/2.,		// half z
									   0.,				// start angle
									   2*M_PI			// total angle
									   );
	eCollimator_log = new G4LogicalVolume(e_collim_tube, Vacuum, "e_collim_log");
	
	// tungsten apertures
	G4Tubs* e_collim_aperture_tube = new G4Tubs("e_collim_aperture_tube", r_inner, r_outer, c_thick/2., 0., 2*M_PI);
	aperture_log = new G4LogicalVolume(e_collim_aperture_tube, Wu, "e_collim_aperture_log");
	for(unsigned int i=0; i<aperture_spacing.size(); i++)
		apertures.push_back(new G4PVPlacement(NULL, G4ThreeVector(0,0,length/2-(aperture_spacing[i]+c_thick/2)),
											  aperture_log, "e_collim_ring_phys_"+itos(i), eCollimator_log, false, i));
}
