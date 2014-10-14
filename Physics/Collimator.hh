#ifndef COLLIMATOR_HH
#define COLLIMATOR_HH

#include <cmath>

/// Base class for calculating aCORN spectrometer acceptance
class EventCollimator {
public:
    /// Constructor
    EventCollimator() { }
    /// Destructor
    virtual ~EventCollimator() { }
    
    /// calculate event pass probability
    virtual double pass() = 0;
    /// transverse magnitude
    static double transv(const double* v) { return sqrt(v[0]*v[0]+v[1]*v[1]); }
    
    double x[3];        ///< vertex position
    double p_e[3];      ///< electron momentum
    double p_p[3];      ///< proton momentum
    
    double pass_e;      ///< probability of electron surviving
    double pass_p;      ///< probability of proton surviving
};

/// Simple collimator (fixed field, hard acceptance)
class SimpleCollimator: public EventCollimator {
public:
    /// Constructor
    SimpleCollimator() { }
    /// Destructor
    virtual ~SimpleCollimator() { }
    
    double B0 = 400;    /// magnetic field [Gauss]
    double r_e = 2.0;   /// electron collimator radius [cm]
    double r_p = 2.5;   /// proton collimator radius [cm]
    
    /// calculate whether event hits wall
    bool hits_wall(double l, const double* p, int sgn);
    
    /// maximum transverse momentum [keV] for collimator of radius r [cm]
    double pt_max(double r) const { return r*B0/3.34; }
                   
    /// calculate event pass probability
    double pass();
};


///////////////////////////////////////////////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

#include "BetaSpectrum.hh"


/// Base class for calculating particle time-of-flight to detector
class ParticleTOF {
public:
    /// Constructor
    ParticleTOF() { }
    /// Destructor
    virtual ~ParticleTOF() { }
    
    /// calculate time-of-flight from initial position [cm] / momentum [keV/c]
    virtual double calcTOF(const double* x, const double* p) const = 0;
};

/// relativistic electron TOF
class ElectronTOF: public ParticleTOF {
public:
    /// Constructor
    ElectronTOF() { }
    
    /// calculate time-of-flight from initial position [cm] / momentum [keV/c]
    virtual double calcTOF(const double* x, const double* p) const;
    
    double det_z = -82;         /// z position of beta detector [cm]
};

/// nonrelativistic proton TOF with electrostatic mirror
class ProtonTOF: public ParticleTOF {
public:
    /// Constructor
    ProtonTOF() { }
    
    /// calculate time-of-flight from initial position [cm] / momentum [keV/c]
    virtual double calcTOF(const double* x, const double* p) const;
    
    double V_mirror = 2.3;      /// potential across electrostatic mirror [kV]
    double L_mirror = 43;       /// length of mirror region [cm]
    double mirror_z = 27.8;     /// start position (0 V) of mirror [cm]
    double det_z = 169.8;       /// z position of proton detector [cm]
};

#endif
