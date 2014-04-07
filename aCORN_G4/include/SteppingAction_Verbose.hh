#ifndef STEPPINGACTION_VERBOSE
#define STEPPINGACTION_VERBOSE

#include <G4SteppingVerbose.hh>

/// class for printing out information at every step
class SteppingAction_Verbose : public G4SteppingVerbose {
public:
	/// constructor
	SteppingAction_Verbose() {}
	/// display info for step
	void StepInfo();
	/// display info at start of track
	void TrackingStarted();
};

#endif
