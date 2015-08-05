/// \file Enums.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef ENUMS_HH
#define ENUMS_HH

#include <utility>
#include <ostream>
#include <Rtypes.h>

#include "AxisEnum.hh"

typedef unsigned int RunNum;
/// wishbone series and run number
typedef std::pair<RunNum, RunNum> RunID;

/// utility for printing std::pair
template<typename T, typename U>
std::ostream& operator<< (std::ostream &out, std::pair<T,U> p) { return out << p.first << "/" << p.second; }

/// simplistic event categorization
enum TriggerCategory {
    TCAT_PROTON = 1<<0, ///< normal proton signal range
    TCAT_PULSER = 1<<1, ///< proton pulser range
    TCAT_TOF_WB = 1<<2  ///< proton TOF covers wishbone range
};

#endif
