/// \file Collimator.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "Collimator.hh"
#include <stdio.h>

/// transverse magnitude
double normTwo(const double* v) { return sqrt(v[0]*v[0]+v[1]*v[1]); }

void LarmorSpiral::setInitial(const double x0[3], const double p[3]) {
    t = 0;
    
    // spiral center offset relative to x0
    cx[0] = -3.34*p[1]/(q*B);
    cx[1] = 3.34*p[0]/(q*B);
    cx[2] = x0[2];
    rL = normTwo(cx);
    phi0 = atan2(-cx[1], -cx[0]);
    gamma = sqrt(1+(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])/(m*m));
    omega = q*B/(gamma*m)*8.987e9;
    
    // convert to absolute position
    for(int i=0; i<3; i++) cx[i] += x0[i];
    rC = normTwo(cx);
}

void LarmorSpiral::calcPos() {
    double phi = phi0 + omega*t;
    xx[0] = cx[0]+rL*cos(phi);
    xx[1] = cx[1]+rL*sin(phi);
    xx[2] = cx[2];
}

void LarmorSpiral::kickMomentum(double dpx, double dpy) {
    // new spiral center
    cx[0] += -3.34*dpy/(q*B);
    cx[1] += 3.34*dpx/(q*B);
    rC = normTwo(cx);
    
    // particle position at deflection relative to center
    double dx[2] = { xx[0]-cx[0], xx[1]-cx[1] };
    
    // new Larmor radius
    rL = normTwo(dx);
    
    // spiral phase CHEAT assuming frequency unchanged (gamma unchanged)
    // TODO recalculate gamma (need to know p_z)
    phi0 = atan2(dx[1], dx[0]) - omega*t;
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////

const double c_light = 2.99792458e10; /// speed of light, cm/s

void ElectronTOF::setInitial(const double x[3], const double p[3]) {
    LarmorSpiral::setInitial(x,p);
    
    vi = p[2]/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+m_e*m_e)/c_light;
}

double ElectronTOF::calcTOF(double z) const {
    if(vi > 0) return INFINITY;
    return (cx[2]-z)*fabs(vi);
}

void ProtonTOF::setInitial(const double x[3], const double p[3]) {
    LarmorSpiral::setInitial(x,p);
    
    double v0 = p[2]/m_p * c_light;                             // initial velocity [cm/s]
    double a = V_mirror/L_mirror/m_p * c_light * c_light;       // acceleration in mirror [cm/s^2]
    t_mr = (-v0 + sqrt(v0*v0 + 2*a*(mirror_z-x[2])))/a;         // ballistic trajectory time in mirror [s]
    
    double E0 = p[2]*p[2]/(2*m_p);                              // initial kinetic energy in z direction [keV/c^2]
    p_exit = sqrt(2*m_p*(E0 + (V_mirror/L_mirror)*(mirror_z-x[2])));   // z momentum exiting mirror [keV/c]
    v_exit = p_exit/m_p * c_light;
}

double ProtonTOF::calcTOF(double z) const {
    double t_det = (z-mirror_z)/v_exit;                         // time from bottom of mirror to detector [s]
    return t_mr + t_det;
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////

double CircleCollimator::pass(ParticleTransport& T) {
    const double rL = T.rL;
    const double rC = T.rC;
    if(rL+rC <= r) return 1;
    if(!n || rL+rC >= r_hard || fabs(rL-rC) >= r) return 0;
    return pow( acos( -(r*r-rL*rL-rC*rC)/(2*rL*rC) ) / M_PI, n );
}
