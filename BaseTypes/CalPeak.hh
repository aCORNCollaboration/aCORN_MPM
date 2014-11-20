#ifndef CALPEAK_HH
#define CALPEAK_HH

#include <string>
using std::string;

#include "FloatErr.hh"

/// calibration peak data
struct CalPeak {
    unsigned int series;
    string pktype;
    string dttype;
    float_err center;
    float_err sigma;
    float_err height;
};

#endif
