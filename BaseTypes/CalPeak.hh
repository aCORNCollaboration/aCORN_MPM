/// \file CalPeak.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

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
