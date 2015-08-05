/// \file Positioner.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef POSITIONER_HH
#define POSITIONER_HH

#include <cmath>
#include <vector>
using std::vector;

/// Extract position information from PMT hit distribution
class Positioner {
public:
    /// Constructor
    Positioner();
    /// Destructor
    virtual ~Positioner() {}
    
    /// calculate center-of-mass position and spread
    void calcPos(const double* L);
    double px[2];               ///< average x,y position
    double sx[2];               ///< sigma x,y position
    /// radius^2
    double r2() const { return px[0]*px[0]+px[1]*px[1]; }
    /// sigma x+y
    double sigma() const { return sqrt(sx[0]*sx[0]+sx[1]*sx[1]); }
    
    /// draw PMT positions on plot
    void drawPMTs(int color = 6, bool drawnum = false) const;
    
    static const unsigned int N = 19;  ///< number of PMTs
    double v[2][2];             ///< lattice vectors in x,y
    int vpos[2][N];             ///< PMT positions in lattice vectors
    double pos[2][N];           ///< PMT positions in x,y
};

#endif
