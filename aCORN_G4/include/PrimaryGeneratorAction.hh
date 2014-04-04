#ifndef PRIMARYGENERATORACTION_HH
#define PRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4VUserEventInformation.hh"
#include "Construct_All.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "ElectronBindingEnergy.hh"
#include "NuclEvtGen.hh"
#include <vector>
#include <Rtypes.h>
#include <TF1.h>

/// User event information for recording primary event weighting
class PrimEvtWeighting: public G4VUserEventInformation {
public:
	/// constructor
	PrimEvtWeighting(double W): w(W) {}
	/// print info
	void Print() const { G4cout << "Primary weighting " << w <<  G4endl; }
	
	double w;	//< event primary weight
};

using namespace std;
class PrimaryGeneratorMessenger;

/// (uncorrected) beta spectrum probability with arbitrary endpoint (for use with TF1)
double genericBetaSpectrum(double* x, double *par);
	
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
	/// constructor
    PrimaryGeneratorAction(DetectorConstruction*);  
	/// destructor
	~PrimaryGeneratorAction();
	
	/// generate primaries for next event
	void GeneratePrimaries(G4Event*);
	/// set gun type
	void SetGunType(G4String val) { gunType = val; }
	/// set default particle type
	void SetParticleType(G4String val) { particleType = val; }
	/// set positioner type
	void SetPositioner(G4String val) { positioner = val; }

	/// set input events tree
	void SetEventFile(G4String val);
	
private:
	G4ParticleGun* particleGun;					//< particle gun primary event thrower
	DetectorConstruction* myDetector;			//< detector geometry
	PrimaryGeneratorMessenger* gunMessenger;	//< messenger of this class
	G4String gunType;							//< event generator gun to use
	G4String particleType;						//< particle type to throw
	G4String positioner;						//< how to position initial events
	
	EventTreeScanner* ETS;						//< Reader for saved input events
	
	// positioning related variables
	G4ThreeVector posOffset;							//< base positioning offset
	void setVertices(std::vector<NucDecayEvent>& evts);	//< set vertex positions for each primary
	
	/// throw a cluster of events
	void throwEvents(const std::vector<NucDecayEvent>& evts, G4Event* anEvent);

	/// print what the particle gun is set up to do
	void displayGunStatus();
	
	long myseed;	//< random seed for event
	/// initialize random seed for event
	void initEventRandomSeed(G4Event* anEvent);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


