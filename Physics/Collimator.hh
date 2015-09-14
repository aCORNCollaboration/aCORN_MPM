/// \file Collimator.hh Simplified particle transport and collimation model for aCORN
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef COLLIMATOR_HH
#define COLLIMATOR_HH

#include <cmath>
#include "PolarizedBetaAsym.hh" // for m_e, m_p constants

/// Particle transport around Larmor spiral in magnetic field
class LarmorSpiral {
public:
    /// Constructor given particle charge and magnetic field [Gauss] given charge and mass [keV/c^2]
    LarmorSpiral(double B0, double z, double mm): B(B0), q(z), m(mm) { }
    /// Destructor
    virtual ~LarmorSpiral() { }
    
    double B;   ///< magnetic field [Gauss]
    double q;   ///< particle charge (electron charge units)
    double m;   ///< particle mass [keV/c^2]
    
    /// initial setup from position [cm], momentum [keV/c]
    virtual void setInitial(const double x0[3], const double p[3]);
    
    /// propagate forward in time
    virtual void addTime(double dt) { t += dt; phi += omega * dt; }
    /// calculate position xx
    void calcPos();
    
    double cx[3];       ///< spiral center position [cm]
    double xx[3];       ///< current position [cm]
    double rL;          ///< Larmor radius [cm]
    double rC;          ///< distance of spiral center from axis [cm]
    double omega;       ///< angular frequency [radians / s]
    double phi;         ///< phase around spiral [radians]
    double t;           ///< accumulated TOF [s]
};


/// Base class for particle transport model, adding z positioning to Larmor progress
class ParticleTransport: public virtual LarmorSpiral {
public:
    /// Constructor
    ParticleTransport(): LarmorSpiral(0,0,0) { }
    /// Destructor
    virtual ~ParticleTransport() { }
    
    /// calculate time [s] required to reach specified z position [cm]
    virtual double calcTOF(double z) const = 0;
};

/// relativistic electron TOF
class ElectronTOF: public ParticleTransport {
public:
    /// Constructor
    ElectronTOF(double B0): LarmorSpiral(B0, -1., m_e) { }
    
    /// set initial conditions
    virtual void setInitial(const double x[3], const double p[3]);
    /// calculate time [s] required to reach specified z position [cm]
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
    ProtonTOF(double B0): LarmorSpiral(B0, 1., m_p) { }
    
    /// set initial conditions
    virtual void setInitial(const double x[3], const double p[3]);
    /// calculate time [s] required to reach specified z position [cm]
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
    virtual double pass(ParticleTransport& T) = 0;
};

/// Circular collimator with probabilistic acceptance
class CircleCollimator: public EventCollimator {
public:
    /// Constructor
    CircleCollimator(double BB, double rr, double rh, unsigned int nn = 1): B(BB), r(rr), r_hard(rh), n(nn) { }
    
    /// calculate event pass probability
    virtual double pass(ParticleTransport& T);
    
    /// maximum transverse momentum [keV] for collimator of radius r [cm]
    double pt_max() const { return r*B/3.34; }
    
    double B;           ///< magnetic field (Gauss)
    double r;           ///< radius (cm)
    double r_hard;      ///< hard-cut radius (cm)
    unsigned int n;     ///< number of apertures (0 for hard cut)
};


#endif
