/// \file PolarizedBetaAsym.cc
#include "PolarizedBetaAsym.hh"
#include <stdio.h>
#include <cmath>
#include <vector>
#include <map>

#ifdef USE_ROOT_MATH
#include <TMath.h>
#define Gamma(x) TMath::Gamma(x)
#else
#define Gamma(x) tgamma(x)
#endif

double Sirlin_g_a2pi(double KE,double KE0,double m) {
    if(KE<=0 || KE>=KE0)
        return 0;
    double b = beta(KE,m);
    double E = KE+m;
    double E0 = KE0+m;
    double athb = atanh(b);
    return (3.*log(m_p/m)-3./4.
    +4.*(athb/b-1.)*((E0-E)/(3.*E)-3./2.+log(2.*(E0-E)/m))
    +4./b*SpenceL(2.*b/(1.+b))
    +athb/b*(2.*(1.+b*b)+(E0-E)*(E0-E)/(6.*E*E)-4.*athb)
    )*alpha/(2.*M_PI);
}

double shann_h_a2pi(double KE, double KE0, double m) {
    if(KE<=0 || KE>=KE0)
        return 0;
    double b = beta(KE,m);
    double E = KE+m;
    double E0 = KE0+m;
    double athb = atanh(b);
    return (3.*log(m_p/m)-3./4.
    + 4.*(athb/b-1.)*((E0-E)/(3.*E*b*b)+(E0-E)*(E0-E)/(24*E*E*b*b)-3./2.+log(2.*(E0-E)/m))
    + 4./b*SpenceL(2.*b/(1.+b))+4*athb/b*(1-athb)
    )*alpha/(2.*M_PI);
}

double Wilkinson_g_a2pi(double W, double W0, double M) {
    if(W>=W0 || W<=1)
        return 0;
    double b = sqrt(W*W-1)/W;
    double athb = atanh(b);
    double g = (3.*log(M)-3./4.
    +4.*(athb/b-1.)*((W0-W)/(3.*W)-3./2.+log(2))
    +4./b*SpenceL(2.*b/(1.+b))
    +athb/b*(2.*(1.+b*b)+(W0-W)*(W0-W)/(6.*W*W)-4.*athb)
    )*alpha/(2.*M_PI)+pow((W0-W),2*alpha/M_PI*(athb/b-1.))-1.;
    return g==g?g:0;
}


double neutronSpectrumCorrectionFactor(double KE) {
    double W = (KE+m_e)/m_e;
    double c = WilkinsonF0(1,W,neutron_R0);             // Fermi function Coulomb
    c *= WilkinsonL0(1,W,neutron_R0);                   // Nonzero charge radius
    c *= CombinedC(1,W,1.,3.,beta_W0,neutron_R0);       // electron/nucleon nonzero size wavefunction convolution
    c *= WilkinsonQ(1,W,beta_W0,proton_M0);             // Coulomb effect on recoil
    c *= (1.+Wilkinson_g_a2pi(W,beta_W0));              // outer radiative corrections
    c *= (1.+Bilenkii59_RWM(W));                        // recoil + WM for free neutron
    return c;
}

double neutronCorrectedBetaSpectrum(double KE) {
    double W = (KE+m_e)/m_e;
    return plainPhaseSpace(W,beta_W0)*neutronSpectrumCorrectionFactor(KE);
}

double Davidson_C1T(double W, double W0, double Z, double R) {
    double p = sqrt(W*W-1);
    double y = alpha*Z*W/p;
    double a2Z2 = alpha*alpha*Z*Z;
    double S0 = sqrt(1-a2Z2);
    double S1 = sqrt(4-a2Z2);
    const double C = pow(Gamma(0.25),2)/sqrt(8*M_PI*M_PI*M_PI); // "Gauss number"
    double sm = 0;
    for(unsigned int n=1; n<10; n++) sm += 1/(n*(n*n+y*y));
    double A = ( (S1+2)/(2*S0+2) * pow(12*Gamma(2.*S0+1.)/Gamma(2.*S1+1.),2) *
    pow(2*p*R,a2Z2/2) * (pow(1-a2Z2/4,2)+y*y) * (1-a2Z2*C/2+a2Z2*y*y*sm/2) );
    
    return (1+S0)*((W0-W)*(W0-W)+A*(W*W-1))/24;
}

double Langer_Cs137_C2T(double W, double W0) {
    const double k = 0.030;
    return (W0-W)*(W0-W)+k*(W*W-1);
}

double Behrens_l2(double W, double W0, double Z, double R) {
    double p = sqrt(W*W-1);
    double y = alpha*Z*W/p;
    double a2Z2 = alpha*alpha*Z*Z;
    double S0 = sqrt(1-a2Z2);
    double S1 = sqrt(4-a2Z2);
    const double C = pow(Gamma(0.25),2)/sqrt(8*M_PI*M_PI*M_PI); // "Gauss number"
    double sm = 0;
    for(unsigned int n=1; n<10; n++) sm += 1/(n*(n*n+y*y));
    return ( (S1+2)/(2*S0+2) * pow(12*Gamma(2.*S0+1.)/Gamma(2.*S1+1.),2) *
    pow(2*p*R,a2Z2/2) * (pow(1-a2Z2/4,2)+y*y) * (1-a2Z2*C/2+a2Z2*y*y*sm/2) );
}

double Behrens_Cs137_C(double W, double W0) {
    double q2 = (W0-W)*(W0-W);
    double p2 = W*W-1;
    double l2 = Behrens_l2(W,W0,56,pow(137,1./3.)*neutron_R0);
    double a1 = 0.000346865*q2 + 0.00331725*l2*p2 - 0.000050327*q2*W + 0.000155636*l2*p2*W + 0.000114834*q2/W;
    double a2 = -0.00427141*q2 - 0.00645269*l2*p2 + 0.000063321*q2*W - 0.000913829*l2*p2*W - 0.000566409*q2/W + 0.0000576232*l2*p2/W;
    double a3 =   0.0131499*q2 + 0.00313793*l2*p2 + 0.00151806*q2*W  + 0.000741562*l2*p2*W - 0.000865957*q2/W - 0.000272219*l2*p2/W;
    double x = 1.07;
    return a1+a2*x+a3*x*x;
}

//-----------------------------------------------------//


double shann_h_minus_g_a2pi(double W, double W0) {
    if(W>=W0 || W<=1)
        return 0;
    double b = sqrt(W*W-1)/W;
    double athb = atanh(b);
    return ( 4.*(athb/b-1.)*(1/(b*b)-1)*(W0-W)/(3*W)*(1+(W0-W)/(8*W))
    +athb/b*(2.-2.*b*b) - (W0-W)*(W0-W)/(6.*W*W) )*alpha/(2.*M_PI);
}

double WilkinsonACorrection(double W) {
    const double W0 = neutronBetaEp/m_e+1.;
    const double mu = 2.792847356-(-1.91304273);        // mu_p - mu_n = 2.792847356(23) - -1.91304273(45) PDG 2010
    const double A_uM = (lambda+mu)/(lambda*(1.-lambda)*(1.+3.*lambda*lambda)*m_p/m_e);
    const double A_1 = lambda*lambda+2.*lambda/3.-1./3.;
    const double A_2 = -lambda*lambda*lambda-3.*lambda*lambda-5.*lambda/3.+1./3.;
    const double A_3 = 2.*lambda*lambda*(1.-lambda);
    return A_uM*(A_1*W0+A_2*W+A_3/W);
}
