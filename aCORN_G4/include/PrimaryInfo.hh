#ifndef PRIMARYINFO_HH
#define PRIMARYINFO_HH

#include "Rtypes.h"
#include "TObject.h"
#include "TObjString.h"

/// ROOT-friendly info on a primary particle
/// we will save a list of these for each event
class PrimaryInfo: public TObject {
public:
	/// constructor
	PrimaryInfo() {}
	
	Float_t vertex[3];	//< primary vertex position [m]
	Float_t p[3];		//< primary momentum [keV]
	
	Float_t KE;			//< primary KE
	Float_t weight;		//< event generator weight
	Long_t  seed;		//< primary random seed
	
	ClassDef(PrimaryInfo, 1);
};

#endif
