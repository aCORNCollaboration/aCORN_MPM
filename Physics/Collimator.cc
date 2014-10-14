#include "Collimator.hh"
#include <stdio.h>

bool SimpleCollimator::hits_wall(double l, const double* p, int sgn) {
    const double d2 = x[0]*x[0] + x[1]*x[1];
    const double l2 = l*l;
    if(d2 > l2) return false;
    const double pt = transv(p);
    double dcosth = sgn*(x[0]*p[1]-x[1]*p[0])/pt;
    double a_max = (l2-d2)/(l+dcosth)/2;
    double a = pt/B0*3.34;  // Larmor radius, cm
    return  a < a_max;
}

double SimpleCollimator::pass() {
    pass_e = (p_e[2] < 0)*hits_wall(r_e, p_e, -1);
    pass_p = hits_wall(r_p, p_p, 1);
    return pass_e * pass_p;
}

////////////////////////////////////////////

const double c_light = 2.99792458e10; /// speed of light, cm/s

double ElectronTOF::calcTOF(const double* x, const double* p) const {
    if(p[2] > 0) return INFINITY;
    return (x[2]-det_z)*fabs(p[2])/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+m_e*m_e)/c_light;
}

double ProtonTOF::calcTOF(const double* x, const double* p) const {
    double v0 = p[2]/m_p * c_light;                             // initial velocity [cm/s]
    double a = V_mirror/L_mirror/m_p * c_light * c_light;       // acceleration in mirror [cm/s^2]
    double t_mr = (-v0 + sqrt(v0*v0 + 2*a*(mirror_z-x[2])))/a;  // ballistic trajectory time in mirror [s]
    
    double E0 = p[2]*p[2]/(2*m_p);                              // initial kinetic energy in z direction [keV/c^2]
    double p_exit = sqrt(2*m_p*(E0 + (V_mirror/L_mirror)*(mirror_z-x[2])));   // z momentum exiting mirror [keV/c]
    double t_det = (det_z-mirror_z)/(p_exit/m_p * c_light);     // time from bottom of mirror to detector [s]
    
    //printf("z=%.2f\tp=%.2f\tmr=%.2f\tdet=%.2f\n", x[2], p[2],t_mr*1e6,t_det*1e6);
    
    return t_mr + t_det;
}
