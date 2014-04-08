#ifndef RUNSETSCANNER_HH
#define RUNSETSCANNER_HH

#include "TChainScanner.hh"
#include "Enums.hh"

#include "QFile.hh"
#include "TagCounter.hh"

#include <map>
#include <cassert>

/// Class for scanning through files specified by run number
class RunSetScanner: public TChainScanner {
public:
	/// constructor
	RunSetScanner(const std::string& treeName);
	
	/// destructor
	virtual ~RunSetScanner();
	
	/// add run to data (return whether successful)
	virtual bool addRun(RunNum rn);
	/// add list of runs to data; return number successfully added
	unsigned int addRuns(const std::vector<RunNum>& rns);
	/// return path to run .root file
	virtual std::string locateRun(RunNum) { assert(false); return ""; }
	
	/// speedload, keeping track of currently loaded run number
	virtual void speedload(unsigned int e);
	/// subclass this for routines when new run is loaded
	virtual void loadNewRun(RunNum) {}
	/// get run number of current event
	virtual RunNum getRun() const { return evtRun; }
	/// check whether this is simulated data
	virtual bool isSimulated() const { return false; }
		
	/// print info about this scanner
	virtual void display();
	/// write run calibrations info to QFile
	void writeCalInfo(QFile& qout, std::string tag);
	
	RunNum evtRun;					//< run number for current event
	double totalTime;				//< combined length of runs in seconds
	TagCounter<RunNum>	runTimes;	//< times for each run loaded
	
protected:
	
	std::vector<RunNum> runlist;			//< list of loaded runs
};


#endif
