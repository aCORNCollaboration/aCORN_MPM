#ifndef ENUMS_HH
#define ENUMS_HH

#include <utility>
#include <ostream>

typedef unsigned int RunNum;
/// wishbone series and run number
typedef std::pair<RunNum, RunNum> RunID;

/// utility for printing std::pair
template<typename T, typename U>
std::ostream& operator<< (std::ostream &out, std::pair<T,U> p) { return out << p.first << "/" << p.second; }

/// axis directions
enum AxisDirection {
    X_DIRECTION = 0,
    Y_DIRECTION = 1,
    Z_DIRECTION = 2,
    T_DIRECTION = 3
};
/// iteration to next axis
inline AxisDirection& operator++(AxisDirection& d) { return d = AxisDirection(d+1); }

/// simplistic event categorization
enum TriggerCategory {
    TCAT_PROTON = 1<<0, ///< normal proton signal range
    TCAT_PULSER = 1<<1, ///< proton pulser range
    TCAT_TOF_WB = 1<<2  ///< proton TOF covers wishbone range
};

#endif
