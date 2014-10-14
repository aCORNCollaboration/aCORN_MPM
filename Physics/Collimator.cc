#include "Collimator.hh"

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
    pass_e = (p_e[2] > 0)*hits_wall(r_e, p_e, -1);
    pass_p = hits_wall(r_p, p_p, 1);
    return pass_e * pass_p;
}
