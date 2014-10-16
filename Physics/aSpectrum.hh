#include "BetaSpectrum.hh"
#include <cstddef>
#include <cassert>
#include <stdint.h>

/// virtual class for providing random number vectors
class Gluck_MC_Rndm_Src {
public:
    /// constructor
    Gluck_MC_Rndm_Src(): u(u0+3) { }
    /// destructor
    virtual ~Gluck_MC_Rndm_Src() { }
    
    /// get next random u[8] for electron, nu, gamma kinematics
    virtual void next() = 0;
    
    double u0[11];      ///< kinematics random array, with optional first 3 position
    double* u;          ///< kinematics array starting at u0[3]
};

/// Implementation of F. Gl\"uck, Computer Physics Communications 101 (1997) 223--231 
class Gluck_beta_MC {
public:
    /// constructor
    Gluck_beta_MC(Gluck_MC_Rndm_Src* R, double M2_F = 1, double M2_GT = 3):
    zeta(M2_F + lambda*lambda*M2_GT), a( (M2_F - lambda*lambda*M2_GT/3.)/zeta ), myR(R) { assert(R); calc_rho(); }
        
    /// Generate weighted event
    double gen_evt_weighted();
    /// Generate un-weighted event by Neumann rejection
    void gen_evt();
    /// Show "efficiency" of MC (5.22)
    void showEffic();
    /// test calculate hard photon decay probability by MC
    void test_calc_P_H();
    /// set random number source
    void SetRandom(Gluck_MC_Rndm_Src* R) { myR = R; }
    
    const double G_F = 1.1663787e-17;   ///< Fermi coupling constant, [/keV^2]
    const double G2_V = G_F*G_F*0.94920; ///< |G_V|^2 = |V_ud G_F g_V|^2
    const double m = m_n;               ///< initial nucleus mass
    const double m_2 = m_e;             ///< mass of charged lepton
    const double Delta = m_n - m_p;     ///< decay energy m - m_f;
    const double C_S = 0.001;           ///< hard photon production cutoff fraction 
    const double zeta;                  ///< spectrum weighting (2.13)
    const double a;                     ///< a_0 (2.13)
        
    double rho_H;       ///< hard decay rate (5.16)
    double rho_0;       ///< 0^th order decay rate (5.18)
    double rho_VS;      ///< soft photon decay rate (5.20)
    double rho_0VS;     ///< soft decay rate = rho_0 + rho_VS;
    double P_H;         ///< probability of hard photon decay (5.21)

    double E_2;         ///< electron total energy
    double p_2;         ///< electron momentum magnitude
    double n_2[3];      ///< electron unit direction
    double beta;        ///< electron velocity v/c (2.10)
    
    double pt2_max = 0; ///< optional limit on maximum electron transverse momentum
    double c_2_min;     ///< optional minimum electron cos theta
    double c_2_wt;      ///< extra weight for c_2 selection
    
    double E0_1;        ///< antineutrino energy in center-of-mass frame (2.10)
    double E_1;         ///< antineutrino energy minus photon (4.9)
    double p_1;         ///< neutrino momentum magnitude
    double n_1[3];      ///< neutrino momentum direction (5.8)
    
    double K;           ///< hard photon energy (5.5)
    double n_gamma[3];  ///< gamma unit direction (5.8)
    double c_gamma;     ///< cos(theta_gamma) (5.6)
    double c_1, c_2;    ///< phase-space cosines
    double phi_gamma;   ///< gamma azimuth relative to electron
    double phi_1, phi_2;///< phase-space azimuths
    
    double p_f[3];      ///< recoil nucleon momentum
    double mag_p_f;     ///< magnitude of recoil momentum
    
    double M_0;         ///< uncorrected decay amplitude (2.11, 2.12)
    double Mtilde;      ///< (3.2)
    double M_VS;        ///< virtual and soft brem amplitude (3.9)
    double M_BR;        ///< hard brem amplitude (4.4)
    
    double evt_w;       ///< calculated event weight for kinematics
    double evt_w0;      ///< uncorrected spectrum weight
    
protected:
    
    Gluck_MC_Rndm_Src* myR;     ///< random number source
    
    double N;           ///< N(beta) (3.3)
    double omega;       ///< hard photon cutoff energy (3.8)
    double np_2[3];     ///< coordinate transform vector (5.10)
    double npp_2[3];    ///< coordinate transform vector (5.10)
    double V_g;         ///< re-weighting function integral (5.15)
    
    double w_avg;       ///< average value of w, (5.23) and (5.16)
    double Wavg_0VS;    ///< average value of W_0VS (5.23)
    double w_max = 0;
    double Wmax_0VS = 0;
    int64_t n_w = 0;
    int64_t n_W_0VS = 0;
    double sum_w = 0;
    double sum_W_0VS = 0;
    
    /// choose initial un-weighted kinematics
    void propose_kinematics();
    /// calculate electron vector and relative coordinates
    void calc_n_2();
    /// calculate unit vector specified relative to n_2
    void vec_rel_n_2(double c, double phi, double* v) const;
    /// calculate beta and N(beta)
    void calc_beta_N();
    /// calculate proton kinematics
    void calc_proton();
    
    /// corrected spectrum probability for virtual and soft brem photons
    double calc_soft();
    /// corrected spectrum probability for hard photon production
    double calc_hard_brem();
    
    /// virtual and soft brem correction (3.10)
    double z_VS() const;
    /// hard brem correction (4.14)
    double z_H() const;
    /// calculate rho_0VS, rho_H non-gamma-emitting rate
    void calc_rho();
};

