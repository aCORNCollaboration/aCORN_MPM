#ifndef ENUMS_HH
#define ENUMS_HH

typedef unsigned int RunNum;

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
