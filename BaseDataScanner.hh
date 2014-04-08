#ifndef BASEDATASCANNER_HH
#define BASEDATASCANNER_HH

#include "Enums.hh"
#include <RTypes.h>
#include "RunSetScanner.hh"

#include "QFile.hh"
#include "TagCounter.hh"

#include <map>
#include "SMExcept.hh"

const unsigned int NCH_MAX = 32;

/// Generic class for processed data TChains
class BaseDataScanner: public RunSetScanner {
public:
	/// constructor
	BaseDataScanner(const std::string& treeName);
	
	/// get info about current event
	virtual Stringmap evtInfo();
	
	Long_t T_p;			// proton time, 10ns units
	Int_t E_p;			// proton energy
	Int_t nP;				// number of protons
	Int_t T_e2p;			// T_p - EPTime
	Int_t E_e;			// Electron energy
	UInt_t nE;				// number of electron PMTs
	Int_t V;				// number of veto PMTs
	UInt_t DetFired;		// detector firing bit map 0-30
	UInt_t DetPiled;		// pileup flags [e- dead zone][det. 0-30]
	Int_t Idx;				// proton index
	Short_t E_PMT[NCH_MAX];	// individual PMT energy
	Char_t T_PMT[NCH_MAX];	// individual PMT time
	Char_t T_PMT_MED;		// median PMT timing offset
	Char_t Max_PMT;			// PMT with largest signal
	Short_t E_Max_PMT;		// energy of max PMT signal
	Int_t flags;			// Event category flags
	
	/// generate event flags
	virtual void makeFlags();
	
	double physicsWeight;		//< event spectrum re-weighting factor
};

#endif
