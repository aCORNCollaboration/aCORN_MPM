/// \file PolarizedBetaAsym.hh calculations for the beta spectrum and various related corrections
#ifndef POLARIZEDBETAASYM_HH
#define POLARIZEDBETAASYM_HH

#include "UnpolarizedBeta.hh"

// references:
// [0] Wilkinson, Analysis of Neutron Beta-Decay, Nucl. Phys. A 377 (1982) 474-504
// [1] Wilkinson, Evaluation of Beta-Decay I,   NIM A 275 (1989) 378-386
// [2] Wilkinson, Evaluation of Beta-Decay II,  NIM A 290 (1990) 509-515
// [3] Wilkinson, Evaluation of Beta-Decay III, NIM A 335 (1995) 305-309
// [4] Wilkinson, Evaluation of Beta-Decay IV,  NIM A 365 (1995) 203-207
// [5] Wilkinson, Evaluation of Beta-Decay V,   NIM A 365 (1995) 497-507

/// uncorrected asymmetry as a function of kinetic energy [keV]
inline double plainAsymmetry(double KE, double costheta=0.5) { return A0_PDG*beta(KE)*costheta; }

/// Shann's h * alpha/2pi radiative correction
double shann_h_a2pi(double KE, double KE0=neutronBetaEp, double m = m_e);
/// (h-g) * alpha/2pi radiative correction to A
double shann_h_minus_g_a2pi(double W, double W0=beta_W0);
/// Wilkinson weak magnetism + g_V*g_A interference + recoil correction to A [0]
double WilkinsonACorrection(double W);

/// combined order-alpha asymmetry corrections
inline double asymmetryCorrectionFactor(double KE) { double W = (KE+m_e)/m_e; return 1.0+WilkinsonACorrection(W)+shann_h_minus_g_a2pi(W); }
/// corrected asymmetry
inline double correctedAsymmetry(double KE, double costheta=0.5) { return plainAsymmetry(KE,costheta)*asymmetryCorrectionFactor(KE); }

#endif
