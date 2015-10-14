#include "Uncorrelated_3Body.hh"
#include <cassert>

void calc_n(double c, double phi, double n[3]) {
    const double s = sqrt(1-c*c);
    n[0] = s*cos(phi);
    n[1] = s*sin(phi);
    n[2] = c;
}

double Uncorrelated_3Body::gen_evt_weighted() {
    assert(myR);
    myR->next(); // random seed
    
    evt_w = 1;
    
    // electron energy, momentum magnitude, velocity
    E_2 = m_2 + (Delta-m_2)*myR->u[0];
    p_2 = sqrt(E_2*E_2 - m_2*m_2);
    beta = sqrt(1-m_2*m_2/E_2/E_2);

    // electron direction, including transverse momentum limiting
    c_2_min = -1;
    if(pt2_max && p_2 > pt2_max)
        c_2_min = sqrt(1.-pt2_max*pt2_max/(p_2*p_2));
    c_2_wt = (1-c_2_min)/2;
    c_2 = c_2_min + (1-c_2_min)*myR->u[1];
    phi_2 = 2*M_PI*myR->u[2];
    calc_n(c_2, phi_2, n_2);
    
    // neutrino energy, direction
    E_1 = Delta - E_2;
    p_1 = E_1; // massless neutrino approximation
    c_1 = 2*myR->u[3] - 1;
    phi_1 = 2*M_PI*myR->u[4];
    calc_n(c_1, phi_1, n_1);
    
    evt_w = evt_w0 = plainPhaseSpace(E_2/m_2);
    
    // proton kinematics
    mag_p_f = 0;
    for(int i=0; i<3; i++) {
        p_f[i] = -n_1[i]*p_1 - n_2[i]*p_2;
        mag_p_f += p_f[i]*p_f[i];
    }
    mag_p_f = sqrt(mag_p_f);
    
    return evt_w;
}