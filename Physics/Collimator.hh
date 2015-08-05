/// \file Collimator.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef COLLIMATOR_HH
#define COLLIMATOR_HH

#include <cmath>
#include "BetaSpectrum.hh"

/// Base class for particle transport model
class ParticleTransport {
public:
    /// Constructor
    ParticleTransport() { }
    /// Destructor
    virtual ~ParticleTransport() { }
    
    /// set initial conditions
    virtual void setVertex(const double* x, const double* p) = 0;
    /// calculate time to z position
    virtual double calcTOF(double z) const = 0;
};

/// relativistic electron TOF
class ElectronTOF: public ParticleTransport {
public:
    /// Constructor
    ElectronTOF() { }
    
    /// set initial conditions
    virtual void setVertex(const double* x, const double* p);
    /// calculate time-of-flight from initial position [cm] / momentum [keV/c]
    virtual double calcTOF(double z) const;
    
    double z0;          ///< initial z position
    double vi;          ///< 1/velocity
};

const double eDet_z = -82; ///< z position of beta detector [cm]
const double pDet_z = 185; ///< z position of proton detector [cm]

/// nonrelativistic proton TOF with electrostatic mirror
class ProtonTOF: public ParticleTransport {
public:
    /// Constructor
    ProtonTOF() { }
    
    /// set initial conditions
    virtual void setVertex(const double* x, const double* p);
    /// calculate time to z position
    virtual double calcTOF(double z) const;
    
    double t_mr;                ///< time spent in mirror
    double p_exit;              ///< momentum at mirror exit
    
    double V_mirror = 3.0;      ///< potential across electrostatic mirror [kV]
    double L_mirror = 43;       ///< length of mirror region [cm]
    double mirror_z = 27.8;     ///< start position (0 V) of mirror [cm]
};



///////////////////////////////////////////////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////



/// Base class for calculating aCORN spectrometer acceptance
class EventCollimator {
public:
    /// Constructor
    EventCollimator() { }
    /// Destructor
    virtual ~EventCollimator() { }
    
    /// calculate event pass probability
    virtual double pass(const double* x, const double* p) = 0;
    
    /// transverse magnitude
    static double normTwo(const double* v) { return sqrt(v[0]*v[0]+v[1]*v[1]); }
};

/// Circular collimator with probabilistic acceptance
class CircleCollimator: public EventCollimator {
public:
    /// Constructor
    CircleCollimator(double BB, double rr, double rh, unsigned int nn = 1): B(BB), r(rr), r_hard(rh), n(nn) { }
    
    /// calculate event pass probability
    virtual double pass(const double* x, const double* p);
    
    /// maximum transverse momentum [keV] for collimator of radius r [cm]
    double pt_max() const { return r*B/3.34; }
    
    double B;           ///< magnetic field (Gauss)
    double r;           ///< radius (cm)
    double r_hard;      ///< hard-cut radius (cm)
    unsigned int n;     ///< number of apertures (0 for hard cut)
};

/*
/// Timing-dependent stack of circular collimators
class StackCollimator: public EventCollimator {
public:
    /// Constructor
    StackCollimator() { }
    
    /// calculate event pass probability
    virtual double pass(const double* x, const double* p);
    
    const ParticleTransport& PT;
};
*/


#endif
