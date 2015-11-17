/// \file Test_Gluck_MC.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

// make Test_Gluck_MC -j4

#include "UnpolarizedNeutronDecay.hh"

#include "OutputManager.hh"
#include "GraphUtils.hh"
#include "ProgressBar.hh"
#include <TStyle.h>
#include <Math/QuasiRandom.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TPad.h>
#include <cassert>

using namespace ROOT::Math;

class RootRandom: public NKine_Rndm_Src {
public:
    RootRandom(size_t n = 8): nrnd(n) { }
    virtual void next() { R.RndmArray(nrnd+3,u0); }
    virtual size_t n_random() const override { return nrnd; }
    
    const size_t nrnd;
    TRandom3 R;
};

class RootQRandom: public NKine_Rndm_Src {
public:
    RootQRandom(size_t n = 8): nrnd(n), QR(nrnd+3) { }
    virtual void next() { QR.Next(u0); }
    virtual size_t n_random() const override { return nrnd; }
    
    const size_t nrnd;
    QuasiRandomSobol QR;
    //QuasiRandomNiederreiter QR;
};

double dot3x(const double* a, const double* b) {
    double s = 0;
    for(int i=0; i<3; i++) s += a[i]*b[i];
    return s;
}

void Gluck93_Table_V() {
    printf(" c\t\t\t\t    r_enu\n");
    for(double c=0.9; c > -0.95; c -= 0.2) {
        printf("%4.1f",c);
        for(double x=0.2; x < 0.95; x += 0.1) printf("\t%5.2f",Gluck93_r_enu(x,c));
        printf("\n");
    }
    printf(" x");
    for(double x=0.2; x < 0.95; x += 0.1) printf("\t %.1f",x);
    printf("\n");
}

void Gluck93_Table_I() {
    printf("\nr_e");
    for(double x = 0.1; x < 0.95; x += 0.1) {
        //printf("\t%.2f",100*Sirlin_g_a2pi(x*neutronBetaEp, neutronBetaEp));
        printf("\t%.2f",100*Wilkinson_g_a2pi(1+x*(beta_W0-1)));
    }
    printf("\n");
}

void G93asym(Gluck93_Distribution& GD, double E2, double& a0, double& aC) {
    GD.calc_Wenu_0Ca(E2,-1);
    double w0m = GD.Wenu_0;
    double wCm = GD.Wenu_0Ca; 
    //double wCm = GD.W_0C * GD.dEfc_dc; // recoil only
    
    GD.calc_Wenu_0Ca(E2,1);
    double w0p = GD.Wenu_0;
    double wCp = GD.Wenu_0Ca;
    //double wCp = GD.W_0C * GD.dEfc_dc; // recoil only
    
    a0 = (w0p-w0m)/(w0p+w0m);
    aC = (wCp-wCm)/(wCp+wCm);
}

void Gluck93_distrib() {
    Gluck93_Distribution GD;
    double a0,aC;
    for(double E2 = m_e+50; E2 < 1300; E2+=50) {
        G93asym(GD,E2,a0,aC);
        
        double cp = B59_rwm_cxn(E2, 1);
        double cm = B59_rwm_cxn(E2, -1);
        double aB59 = (cp-cm)/(cp+cm);
        
        printf("%g\t%.3f\t%.3f\t%.3f\t\t%.3f\n", E2-m_e, 100*a0, 100*aC, 100*(aC-a0)/a0, 100*aB59/a0);
    }
}

void GluckBetaCompare() {
    OutputManager OM("Gluck_MC_Test", "/home/mpmendenhall/Desktop/aCORN_MC/");
    gStyle->SetOptStat("");
    
    TH1D* hBeta_A[2];
    TH1D* hBeta_MC[2];
    for(int i=0; i<2; i++) {
        hBeta_A[i] = new TH1D(i?"hBeta_A1":"hBeta_A0", "Analytical beta spectrum", 1000, 0, .800);
        hBeta_A[i]->GetXaxis()->SetTitle("energy [MeV]");
        hBeta_A[i]->GetYaxis()->SetTitle("rate [1/decay/MeV]");
        
        hBeta_MC[i] = new TH1D(i?"hBeta_MC1":"hBeta_MC0", "MC beta spectrum", 100, 0, .800);
        hBeta_MC[i]->GetXaxis()->SetTitle("energy [MeV]");
        hBeta_MC[i]->GetYaxis()->SetTitle("rate [1/decay/MeV]");
    }
    
    // analytical calculation
    double sp = 0;
    const int nthpts = 100;
    BetaSpectrumGenerator BSG(1,1,neutronBetaEp);
    for(int i=1; i<=hBeta_A[0]->GetNbinsX(); i++) {
        double KE = 1000*hBeta_A[0]->GetBinCenter(i);
        double p = BSG.decayProb(KE);
        if(!p) continue;
        for(int j=0; j < nthpts; j++) {
            double cos_thn = -1 + 2.*(j+0.5)/nthpts;
            double awt = GM78_radiative_cxn(KE+m_e, cos_thn)*B59_rwm_cxn(KE + m_e, cos_thn)/(1.+Bilenkii59_RWM(KE/m_e+1));
            hBeta_A[cos_thn > 0]->Fill(KE/1000., p*awt);
            sp += p*awt;
        }
    }
    for(int i=0; i<2; i++) {
        hBeta_A[i]->Rebin(10);
        normalize_to_bin_width(hBeta_A[i], 1./sp);
    }
    
    RootQRandom RQR;
    Gluck_beta_MC G(&RQR);
    //G.test_calc_P_H();
    
    /// weighted simulation mode
    int nToSim = 1e8;
    ProgressBar* PB = new ProgressBar(nToSim, 50);
    double sw = 0;
    for(int i=0; i<nToSim; i++) {
        PB->update(i);
        G.gen_evt_weighted();
        if(G.evt_w <= 0) continue;
        
        G.evt_w *= G.coulomb_cxn()*G.rwm_cxn();
        sw += G.evt_w;
        
        double cos_thn = dot3x(G.n_1, G.n_2);
        hBeta_MC[cos_thn > 0]->Fill((G.E_2 - m_e)/1000, G.evt_w);
    }
    delete PB;
    printf(" Done.\n");
    G.showEffic();
    for(int i=0; i<2; i++) {
        normalize_to_bin_width(hBeta_MC[i], 1./sw);
    }
    
    for(int i=0; i<2; i++) hBeta_A[i]->SetLineColor(2);
    hBeta_A[0]->Draw("HIST");
    hBeta_A[1]->Draw("HIST Same");
    for(int i=0; i<2; i++) hBeta_MC[i]->Draw("HIST Same");
    OM.printCanvas("hBeta");
    
    hBeta_A[0]->Add(hBeta_A[1], -1);
    hBeta_MC[0]->Add(hBeta_MC[1], -1);
    hBeta_MC[0]->Add(hBeta_A[0], -1.);
    hBeta_MC[0]->Draw("HIST");
    OM.printCanvas("hBeta_Resid");
    
    printf("Integral error: %g\n", hBeta_MC[0]->Integral("width"));
}

void GluckSirlinCompare() {
    RootQRandom RQR;
    Gluck_beta_MC G(&RQR);
    for(double E2 = 550; E2 < 782+m_e; E2 += 50) {
        double gS = Sirlin_g_a2pi(E2-m_e, neutronBetaEp);
        double gW = Wilkinson_g_a2pi(E2/m_e);
        double gG = G.recalc_Sirlin_g_a2pi(E2);
        printf("%g\t%g\t%g\t%g\n", E2, gS, gW, gG);
    }
}


int main(int, char**) {
 
    //Gluck93_Table_V();
    //Gluck93_Table_I();
    //Gluck93_distrib();
    //GluckBetaCompare();
    GluckSirlinCompare();
    
    
    return EXIT_SUCCESS;
    
}
