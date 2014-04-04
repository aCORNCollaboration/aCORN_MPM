#ifndef PHYSICSLIST_LIVERMORE_hh
#define PHYSICSLIST_LIVERMORE_hh

#include "G4VUserPhysicsList.hh"

/// UCNA electromagnetic physics list based on Livermore low-energy routines
class PhysicsList_Livermore: public G4VUserPhysicsList {
public:
	/// constructor
    PhysicsList_Livermore() {}
	/// set particle cuts
	virtual void SetCuts();
	
protected:
    /// list particles to consider
    virtual void ConstructParticle();
	/// list physics processes to consider
    virtual void ConstructProcess();
    /// electromagnetic physics processes
	virtual void StandardEM();
};

#endif
