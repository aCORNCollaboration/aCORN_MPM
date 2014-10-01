#include "aSpectrum.hh"
#include <stdio.h>

double dot3(const double* a, const double* b) {
    double s = 0;
    for(int i=0; i<3; i++) s += a[i]*b[i];
    return s;
}

double Gluck_beta_MC::z_VS() const {
    // (3.10)
    if(!omega) return 0;
    return alpha/M_PI * (3./2.*log(m_p/m_2) + 2*(N/beta-1)*log(2*omega/m_2)
                         +2*N/beta*(1-N) + 2/beta*SpenceL(2*beta/(1+beta)) - 3./8.);
}

double Gluck_beta_MC::z_H() const {
    // (4.14), noting E0_1/omega = 1/C_S
    return alpha/M_PI * (2*(N/beta-1)*(log(1/C_S) + E0_1/(3*E_2) -3./2.)
                         + N/beta * E0_1*E0_1/(12*E_2*E_2));
    
}

double Gluck_beta_MC::calc_soft() {
    
    vec_rel_n_2(c_1, phi_1, n_1);
  
    // (2.12) 
    M_0 = 16 * G2_V * zeta * m*m * E0_1 * E_2 * (1 + a * beta * c_1); 
    
    // (3.2)
    Mtilde = -alpha/M_PI * 16 * G2_V * (1-beta*beta)/beta * N * m*m * E0_1 * E_2 * zeta;
    
    // (3.9)
    M_VS = z_VS() * M_0 + Mtilde;
    
    // weight function
    double W_0VS = beta * E0_1 * E_2 * (M_0 + M_VS);
    evt_w0 =  beta * E0_1 * E_2 * M_0 * evt_w;
    if(W_0VS > Wmax_0VS) Wmax_0VS = W_0VS; 
    sum_W_0VS += W_0VS;
    n_W_0VS++;
    return W_0VS;
}

void Gluck_beta_MC::vec_rel_n_2(double c, double phi, double* v) const {
    // (5.11)
    const double s = sqrt(1-c*c);
    // (5.9)
    double n_perp[3];
    for(int i=0; i<3; i++) n_perp[i] = np_2[i]*cos(phi) + npp_2[i]*sin(phi);
    // (5.8)
    for(int i=0; i<3; i++) v[i] = n_2[i]*c + n_perp[i] * s;
}
    
double Gluck_beta_MC::calc_hard_brem() {
    
    // (5.11)
    const double s_1 = sqrt(1-c_1*c_1);
    // (5.12)
    n_1[0] = s_1*cos(phi_1); n_1[1] = s_1*sin(phi_1); n_1[2] = c_1;
    
    // (5.5)
    K = omega * exp(-myR->u[5]*log(C_S));
    
    // (4.9)
    E_1 = Delta - E_2 - K;
    
    // (5.6)
    c_gamma = (1-(1+beta)*exp(-2*N*myR->u[6]))/beta;
    
    // (5.7)
    phi_gamma = 2*M_PI*myR->u[7];
    
    vec_rel_n_2(c_gamma, phi_gamma, n_gamma);
        
    // (5.13)
    const double n_1_dot_n_2 = dot3(n_1, n_2);
    const double n_1_dot_n_gamma = dot3(n_1, n_gamma);
    const double p_1_dot_k = E_1*K*n_1_dot_n_gamma;
    const double p_2_dot_k = beta*E_2*K*c_gamma;
    const double p_1_dot_p_2 = beta*E_1*E_2*n_1_dot_n_2;
    // (4.8)
    const double p4_2_dot_k4 = E_2*K - p_2_dot_k; 
    
    // (5.3) approximate distribution function
    const double g = beta*E_2/(2*N*p_2_dot_k);
    //printf("\tc_gamma = %g g = %g\n",c_gamma,g);
    
    // (4.7)
    const double Psq = 1./K/K + m_2*m_2/p4_2_dot_k4/p4_2_dot_k4 - 2*E_2/K/p4_2_dot_k4;
    // (4.5)
    const double H_0 = E_1*(-(E_2+K)*Psq + K/p4_2_dot_k4);
    // (4.6)
    const double H_1 = ( p_1_dot_p_2 * ( -Psq + 1./p4_2_dot_k4 )
                        + p_1_dot_k * ( (E_2+K)/K -m_2*m_2/p4_2_dot_k4)/p4_2_dot_k4);
    // (4.3)
    const double esq = 4*M_PI*alpha;
    // (4.4)
    M_BR = 16*G2_V*zeta*m*m*esq*(H_0 + a*H_1);
    
    // (5.14)
    double w = (K * beta * E_1 * E_2 * M_BR)/(pow(2,13)*pow(M_PI,8)*m*m*g);
    if(w > w_max) w_max = w;
    sum_w += w;
    n_w++;
    return w;
}

void Gluck_beta_MC::calc_beta_N() {
    // (2.10)
    beta = sqrt(1-m_2*m_2/E_2/E_2);
    if(!(beta==beta)) beta = 0;
    // (3.3)
    N = 0.5*log((1+beta)/(1-beta));
    // (3.8) soft brem cutoff
    omega = C_S * E0_1;
}

void Gluck_beta_MC::propose_kinematics() {
    
    K = 0;
    evt_w = 1;
    
    // (5.4)
    E_2 = m_2 + (Delta-m_2)*myR->u[0];
    p_2 = sqrt(E_2*E_2 - m_2*m_2);
    
    c_2_min = -1;
    if(pt2_max && p_2 > pt2_max)
        c_2_min = sqrt(1.-pt2_max*pt2_max/(p_2*p_2));
    c_2_wt = (1-c_2_min)/2;
    c_2 = c_2_min + (1-c_2_min)*myR->u[2];

    // (2.10)
    E_1 = E0_1 = Delta - E_2;
       
    // (5.7)
    c_1 = 2*myR->u[1] - 1;
    phi_1 = 2*M_PI*myR->u[3];
    phi_2 = 2*M_PI*myR->u[4];
     
    calc_beta_N();
    calc_n_2();
}

void Gluck_beta_MC::calc_n_2() {
    // (5.11)
    const double s_2 = sqrt(1-c_2*c_2);
    // (5.12)
    n_2[0] = s_2*cos(phi_2); n_2[1] = s_2*sin(phi_2); n_2[2] = c_2;
    // (5.10)
    np_2[0] = -sin(phi_2);      np_2[1] = cos(phi_2);           np_2[2] = 0;
    npp_2[0] = -c_2*cos(phi_2); npp_2[1] = -c_2*sin(phi_2);     npp_2[2] = s_2;
}

double Gluck_beta_MC::gen_evt_weighted() {
    assert(myR);

    if(P_H < myR->selectBranch()) {
        myR->next_0();
        propose_kinematics();
        evt_w *= calc_soft()/Wavg_0VS;
        evt_w0 /= Wavg_0VS;
    } else {
        myR->next_H();
        propose_kinematics();
        evt_w *= calc_hard_brem()/w_avg;
        evt_w0 = 0;
    }
    calc_proton();
    
    return evt_w;
}

void Gluck_beta_MC::gen_evt() {
    assert(myR);
    evt_w = 1;
    
    if(P_H < myR->selectBranch()) {
        do {
            myR->next_0();
            propose_kinematics();
        } while(calc_soft()/Wmax_0VS < myR->selectBranch());
    } else {
        do {
            myR->next_H();
            propose_kinematics();
        } while(calc_hard_brem()/w_max < myR->selectBranch());
    }
    calc_proton();
}

void Gluck_beta_MC::calc_proton() {
    p_1 = E_1; // massless neutrino approximation
    for(int i=0; i<3; i++) p_f[i] = -n_1[i]*p_1 - n_2[i]*p_2 - K*n_gamma[i];
}

void Gluck_beta_MC::calc_rho() {
    
    printf("Calculating correction rates...\n");
    
    // numerical integration by Simpson's Rule
    const int npts = 1001;
    rho_0 = rho_VS = rho_H = 0;
    const double C = G2_V * zeta / (2*pow(M_PI,3));
    for(int i=0; i<=npts; i++) {
        
        E_2 = m_2 + i*(Delta-m_2)/npts;
        E0_1 = Delta - E_2;
        calc_beta_N();
        
        if(!beta) continue;
        //printf("E_2 = %g, beta = %g, z_VS = %g, z_H = %g\n", E_2, beta, z_VS(), z_H());
        
        const int scoeff = (i==0 || i==npts)? 1: (i%2)? 4:2;
        // (5.19)
        const double w0 = C * beta * E0_1*E0_1 * E_2*E_2;
        // (5.18)
        rho_0 += scoeff * w0;
        // (5.20)
        rho_VS += scoeff * w0 * (z_VS() - alpha/M_PI * N * (1-beta*beta)/beta);
        // (4.13)
        rho_H += scoeff * C * beta * E0_1*E0_1 * E_2*E_2 * z_H(); 
    }
    const double nrm = (Delta-m_2)/npts/3.;
    rho_0 *= nrm;
    rho_VS *= nrm;
    rho_H *= nrm;
    
    rho_0VS = rho_0 + rho_VS;
    
    // (5.21)
    P_H = rho_H/(rho_H + rho_0VS);
    
    // (5.15)
    V_g = -32*pow(M_PI,3)*(Delta-m_2)*log(C_S);
    
    w_avg = rho_H / V_g;
    
    printf("\trho_0 = %g, rho_VS = %g, rho_H = %g => P_H = %g\n", rho_0, rho_VS, rho_H, P_H);
}

void Gluck_beta_MC::test_calc_P_H() {
    
    // (5.16), (5.17)
    double t_rho_H = 0;
    double sw2 = 0;
    const int n_H = 1000000;
    printf("Calculating P_H using %g points... ",double(n_H));
    for(int i=0; i<n_H; i++) {
        if(!(i%(n_H/20))) { printf("*"); fflush(stdout); }
        
        myR->next_H();
        propose_kinematics();
        double w = calc_hard_brem();
        t_rho_H +=  w;
        sw2 += w*w;
        
        myR->next_0();
        propose_kinematics();
        calc_soft();
    }
    double dt_rho_H = sqrt(sw2 - t_rho_H*t_rho_H/n_H)*V_g/n_H;
    t_rho_H *= V_g/n_H;
    
    // (5.21)
    double t_P_H = t_rho_H/(rho_H + rho_0VS);
    printf("\tMC rho_H = %g +/- %g, P_H = %g\n", t_rho_H, dt_rho_H, t_P_H);
    
    showEffic();
}

void Gluck_beta_MC::showEffic() {
    // (5.22)
    if(n_W_0VS) { Wavg_0VS = sum_W_0VS/n_W_0VS; printf("Wmax_0VS = %g; E_0VS = %g %%\n", Wmax_0VS, 100*sum_W_0VS/n_W_0VS/Wmax_0VS); }
    if(n_w) { w_avg = sum_w/n_w; printf("w_max = %g; E_H = %g %%\n", w_max, 100*w_avg/w_max); }
}
