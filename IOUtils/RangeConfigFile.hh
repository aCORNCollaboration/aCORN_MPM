/// \file RangeConfigFile.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef RANGECONFIGFILE_HH
#define RANGECONFIGFILE_HH

#include "SMFile.hh"

/// class for looking up manually-stored miscellaneous analysis info
class RangeConfigFile: public SMFile {
public:
    /// constructor
    RangeConfigFile(string fname): SMFile(fname) {}
    
    /// get (double start,double end) pairs for key (or any other named pairs)
    vector< std::pair<double,double> > getRanges(const std::string& key, const std::string& k1="start", const std::string& k2="end") const;
    
    /// get matching keys for item in range
    vector<Stringmap> getInRange(const std::string& key,
                                      const double x,
                                      const std::string& k1="runStart",
                                      const std::string& k2="runEnd") const;
    static RangeConfigFile MI;       ///< static global instance to use
};

/// simple class for cuts from Stringmap
class RangeCut {
public:
    /// constructor
    RangeCut(const Stringmap& m = Stringmap());
    /// constructor with start, end times
    RangeCut(double s, double e): start(s), end(e) {}
    
    /// check if value is in range
    inline bool inRange(double x) const { return start <= x && x <= end; }
    
    double start;       ///< cut minimum
    double end;         ///< cut maximum
};

/// simple class for value + cuts range
class CutVariable {
public:
    /// constructor
    CutVariable(string sn=""): sname(sn) {}
    /// check if in range
    inline bool inRange() const { return R.inRange(val); }
    string sname;  ///< sensor name (for pedestal subtraction)
    double val;         ///< stored value
    RangeCut R;         ///< cuts range
};

Stringmap loadCut(RunNum rn, const std::string& cutName);

void loadRangeCut(RunNum rn, CutVariable& c, const std::string& cutName);

#endif
