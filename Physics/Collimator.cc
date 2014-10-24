#include "Collimator.hh"
#include <stdio.h>

double CircleCollimator::pass(const double* x, const double* p) {
    const double aa[2] = { -3.34*p[1]/B, 3.34*p[0]/B }; // offset to spiral center [cm]
    const double a = normTwo(aa);                       // Larmor radius [cm]
    const double cc[2] = { x[0] + aa[0], x[1] + aa[1] };// position of spiral center [cm]
    const double c = normTwo(cc);                       // center distance from origin [cm]
    if(a+c <= r) return 1;
    if(!n || a+c >= r_hard || fabs(a-c) >= r) return 0;
    return pow( acos( -(r*r-a*a-c*c)/(2*a*c) ) / M_PI, n );
}

////////////////////////////////////////////

const double c_light = 2.99792458e10; /// speed of light, cm/s

void ElectronTOF::setVertex(const double* x, const double* p) {
    z0 = x[2];
    vi = p[2]/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+m_e*m_e)/c_light;
}

double ElectronTOF::calcTOF(double z) const {
    if(vi > 0) return INFINITY;
    return (z0-z)*fabs(vi);
}

void ProtonTOF::setVertex(const double* x, const double* p) {
    double v0 = p[2]/m_p * c_light;                             // initial velocity [cm/s]
    double a = V_mirror/L_mirror/m_p * c_light * c_light;       // acceleration in mirror [cm/s^2]
    t_mr = (-v0 + sqrt(v0*v0 + 2*a*(mirror_z-x[2])))/a;         // ballistic trajectory time in mirror [s]
    
    double E0 = p[2]*p[2]/(2*m_p);                              // initial kinetic energy in z direction [keV/c^2]
    p_exit = sqrt(2*m_p*(E0 + (V_mirror/L_mirror)*(mirror_z-x[2])));   // z momentum exiting mirror [keV/c]
}

double ProtonTOF::calcTOF(double z) const {
    double t_det = (z-mirror_z)/(p_exit/m_p * c_light);     // time from bottom of mirror to detector [s]
    return t_mr + t_det;
}
