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
    virtual double pass(const double* x, const double* p) = 0;
    
    /// transverse magnitude
    static double normTwo(const double* v) { return sqrt(v[0]*v[0]+v[1]*v[1]); }
};

/// Simple circular collimator (fixed field, hard acceptance)
class HardCollimator: public EventCollimator {
public:
    /// Constructor
    HardCollimator(double BB, double rr): B(BB), r(rr) { }
    /// Destructor
    virtual ~HardCollimator() { }
    
    /// calculate event pass probability
    virtual double pass(const double* x, const double* p);
        
    /// maximum transverse momentum [keV] for collimator of radius r [cm]
    double pt_max() const { return r*B/3.34; }
    
    double B;   ///< magnetic field (Gauss)
    double r;   ///< radius (cm)
};

/// Circular collimator with probabilistic acceptance
class SoftCollimator: public HardCollimator {
public:
    /// Constructor
    SoftCollimator(double BB, double rr, unsigned int nn = 1): HardCollimator(BB,rr), n(nn) { }
    
    /// calculate event pass probability
    virtual double pass(const double* x, const double* p);
    
    unsigned int n;     ///< number of apertures
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
    
    double V_mirror = 3.0;      /// potential across electrostatic mirror [kV]
    double L_mirror = 43;       /// length of mirror region [cm]
    double mirror_z = 27.8;     /// start position (0 V) of mirror [cm]
    double det_z = 169.8;       /// z position of proton detector [cm]
};

#endif
