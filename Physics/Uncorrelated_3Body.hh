/// \file Uncorrelated_3Body.hh Simplistic test neutron beta decay generator with uncorrelated electron, neutrino momenta
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef UNCORRELATED_3BODY_HH
#define UNCORRELATED_3BODY_HH

#include "UnpolarizedBeta.hh"

/// virtual class for providing random number vectors
class U3Body_Rndm_Src {
public:
    /// constructor
    U3Body_Rndm_Src(): u(u0+3) { }
    /// destructor
    virtual ~U3Body_Rndm_Src() { }
    
    /// get next random u[8] for electron, nu, gamma kinematics
    virtual void next() = 0;
    
    double u0[8];       ///< kinematics random array, with optional first 3 position
    double* u;          ///< kinematics array starting at u0[3]
};

class Uncorrelated_3Body {
public:
    /// constructor
    Uncorrelated_3Body(U3Body_Rndm_Src* R): myR(R) { }
        
    /// Generate weighted event
    double gen_evt_weighted();
    
    /// calculate cos theta between electron, neutrino
    double cos_theta_e_nu() const { return n_1[0]*n_2[0] + n_1[1]*n_2[1] + n_1[2]*n_2[2]; }
    
    void showEffic() const { }
    
    const double m = m_n;               ///< initial nucleus mass
    const double m_2 = m_e;             ///< mass of charged lepton
    const double Delta = m_n - m_p;     ///< decay energy m - m_f;

    double E_2;         ///< electron total energy [keV]
    double p_2;         ///< electron momentum magnitude [keV/c]
    double n_2[3];      ///< electron unit direction
    double beta;        ///< electron velocity v/c (2.10)
    
    double pt2_max = 0; ///< optional limit on maximum electron transverse momentum
    double c_2_min;     ///< optional minimum electron cos theta
    double c_2_wt;      ///< extra weight for c_2 selection
    
    double E0_1;        ///< antineutrino energy in center-of-mass frame [keV] (2.10)
    double E_1;         ///< antineutrino energy minus photon [keV] (4.9)
    double p_1;         ///< neutrino momentum magnitude [keV/c]
    double n_1[3];      ///< neutrino momentum unit direction (5.8)
    
    double c_1, c_2;    ///< phase-space cosines
    double phi_gamma;   ///< gamma azimuth relative to electron
    double phi_1, phi_2;///< phase-space azimuths
    
    double p_f[3];      ///< recoil nucleon momentum
    double mag_p_f;     ///< magnitude of recoil momentum
    
    double evt_w;       ///< calculated event weight for kinematics
    double evt_w0;      ///< weight coming from uncorrected spectrum shape   
    
    U3Body_Rndm_Src* myR;       ///< random number source
    
    double K = 0;       ///< unused radiative decay member
};

#endif
