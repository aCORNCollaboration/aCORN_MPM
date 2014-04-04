#ifndef PHYSICSLIST_495
#define PHYSICSLIST_495

#include "G4VModularPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"

class PhysicsList_495: public G4VModularPhysicsList {
public:
	
	PhysicsList_495(bool usePenelope);
	virtual ~PhysicsList_495();
	
	void ConstructParticle();
    
	void SetCuts();
	void SetCutForGamma(G4double);
	void SetCutForElectron(G4double);
	void SetCutForPositron(G4double);
	
	void ConstructProcess();
	
private:
		
	G4double cutForGamma;
	G4double cutForElectron;
	G4double cutForPositron;
    
	G4String emName;
	G4VPhysicsConstructor* emPhysicsList;
};

#endif
