/// \file SpectrumShape.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "UnpolarizedNeutronDecay.hh"
#include "Collimator.hh"
#include "GraphicsUtils.hh"
#include "StringManip.hh"
#include "PathUtils.hh"
#include "OutputManager.hh"
#include "ProgressBar.hh"

#include <TStyle.h>
#include <Math/QuasiRandom.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TPad.h>
#include <cassert>
#include <TF1.h>

#include <gsl/gsl_qrng.h>

using namespace ROOT::Math;

class RootRandom: public NKine_Rndm_Src {
public:
    RootRandom(size_t n = 8): nrnd(n) { }
    virtual void next() override { R.RndmArray(nrnd+3,u0); }
    virtual size_t n_random() const override { return nrnd; }
    
    const size_t nrnd;
    TRandom3 R;
};

class RootQRandom: public NKine_Rndm_Src {
public:
    RootQRandom(size_t n = 8): nrnd(n), QR(nrnd+3) { }
    virtual void next() override { QR.Next(u0); }
    virtual size_t n_random() const override { return nrnd; }
    
    const size_t nrnd;
    QuasiRandomSobol QR;
    //QuasiRandomNiederreiter QR;
};

class GSLQRandom: public NKine_Rndm_Src {
public:
    GSLQRandom(size_t n = 8): nrnd(n) {
        // choices:
        // gsl_qrng_niederreiter_2: fails for n+3 = 11
        // gsl_qrng_sobol
        // gsl_qrng_halton
        // gsl_qrng_reversehalton 
        myQRNG = gsl_qrng_alloc(gsl_qrng_sobol, n+3);
    }
    ~GSLQRandom() { gsl_qrng_free(myQRNG); }
    virtual void next() override { gsl_qrng_get(myQRNG, u0); }
    virtual size_t n_random() const override { return nrnd; }
    
    const size_t nrnd;
    gsl_qrng* myQRNG;
};

/// 3-body decay generator that pairs "flipped" neutrino directions
class NuFlipper: public N3BodyUncorrelated {
public:
    /// Constructor
    NuFlipper(NKine_Rndm_Src* R): N3BodyUncorrelated(R) { }
    
    /// Generate weighted event
    virtual void gen_evt_weighted() override {
        if(!(nflip++ % 2)) N3BodyUncorrelated::gen_evt_weighted();
        else {
            n_1[2] *= -1;
            calc_proton();
        }
    }
    
    size_t nflip = 0;
};

void circle_the_square(double* x, double r0) {
    double phi = x[0]*2*M_PI;
    double r = sqrt(x[1]*r0*r0);
    x[0] = r*cos(phi);
    x[1] = r*sin(phi);
}

TH1* calcAsym(TH1* h1, TH1* h2) {
    TH1* hAsym = (TH1*)h1->Clone("asym");
    for(int i=0; i <= h1->GetNbinsX()+1; i++) {
        double n1 = h1->GetBinContent(i);
        double n2 = h2->GetBinContent(i);
        if(n1+n2)
            hAsym->SetBinContent(i,(n1-n2)/(n1+n2));
    }
    return hAsym;
}

/// Calculations for electric field momentum deflection
class MomentumDeflector {
public:
    /// Constructor
    MomentumDeflector() { }
    
    /// momentum impulse in x direction due to wires [cm V/cm]
    double wires_momentum_deflection(double x) const {
        return -(floor(x/d) + 0.5 - x/d)*d*E0*0.5;
    }
    
    /// momentum deflection from endcap mismatch [cm V/cm]
    void radial_momentum_deflection(double x, double y, double& Vx, double& Vy) const {
        double rr = x*x + y*y;
        double r = sqrt(rr);
        double rrrr = rr*rr;
        double V = 1.1703106983790565*r -0.54657810070685897*rr +0.51766136215597991*r*rr -0.18217095995451557*rrrr +0.025820766945784723*r*rrrr;
        Vx = x*V/r;
        Vy = y*V/r;
    }
    
    double E0 = 70;     ///< mirror field [V/cm]
    double d = 0.19;    ///< wire spacing
};

int main(int, char**) {
    setupSlideStyle();
    
    OutputManager OM("Simulated", getEnvSafe("ACORN_SUMMARY")+"/Simulated_test/");
    
    // kinematics generator
    GSLQRandom RQR(5);
    
    //Gluck_beta_MC G(&RQR);
    //G.test_calc_P_H();
    //G.P_H.q = 0.5;
    
    N3BodyUncorrelated G(&RQR);
    
    const double a0 = calc_a0();
    
    assert(G.n_random() <= RQR.n_random());
    
    const double B0 = 364;      // field [Gauss]
    const double r_beam = 3.0;  // beam radius [cm]
    ElectronTOF eTOF(B0);
    ProtonTOF pTOF(B0);
    CircleCollimator eCol(B0, 2.75, 5.08, 0);
    CircleCollimator pCol(B0, 4.0, 5.08, 0);
    MomentumDeflector MD;
    
    G.pt2_max = eCol.pt_max();  // limit max electron transverse momentum to save calculation time
    
    ////////////// Wishbone timing estimates
    const double zmin[3] = {0,0,-r_beam};
    const double zmax[3] = {0,0,r_beam};
    const double pemax[3] = {0, 0, sqrt(G.Delta*G.Delta - m_e*m_e)};
    const double pnumax[3] = {0, 0, -(G.Delta-m_e)};
    pTOF.setInitial(zmin, pnumax);
    const double tof_max = pTOF.calcTOF(pDet_z);
    pTOF.setInitial(zmax, pemax);
    const double tof_min = pTOF.calcTOF(pDet_z);
    double tof_mid = 0.5*(tof_max+tof_min);
    printf("Min proton TOF %.2f us, mid %.2f us, max %.2f us\n", tof_min*1e6, tof_mid*1e6, tof_max*1e6);
    tof_mid = 3.4e-6;

    TH1D* hSpec = OM.registeredTH1D("hSpec","Corrected beta spectrum",200,0,800);
    
    TH1D* haCorn[2][2]; // energy spectrum for [branch][radiative corrected?]
    for(int i=0; i<2; i++) {
        for(int j=0; j<2; j++) {
            haCorn[i][j] = OM.registeredTH1D("haCorn_"+to_str(i)+"_"+to_str(j), "asymmetry beta spectrum", 200,0,800);
            haCorn[i][j]->SetLineColor(2+2*i);
        }
    }
    
    TH2F* hPassPos = OM.registeredTH2F("hPassPos","Vertex passed collimator",50,-4,4,50,-4,4);
    
    TH2F* hWishbone[2]; // wishbone for upper/lower branches
    for(int i=0; i<2; i++) {
        hWishbone[i] = OM.registeredTH2F("hWishbone_"+to_str(i),"simulated wishbone", 400, 0, 800, 400, 2, 5);
        hWishbone[i]->GetXaxis()->SetTitle("Electron energy [keV]");
        hWishbone[i]->GetYaxis()->SetTitle("Proton TOF [#mus]");
    }
    
    TH2F* hAccept = OM.registeredTH2F("hAccept","simulated acceptance", 200, 0, 800, 200, -1, 1);
    hAccept->GetXaxis()->SetTitle("Electron energy [keV]");
    hAccept->GetYaxis()->SetTitle("cos #theta_{e#nu}");
    
    TH2F* hpAccept = OM.registeredTH2F("hpAccept","simulated proton acceptance", 200, 0, 800, 200, -1, 1);
    hpAccept->GetXaxis()->SetTitle("Electron energy [keV]");
    hpAccept->GetYaxis()->SetTitle("cos #theta_{p}");
    
    TH2F* hnuAccept = OM.registeredTH2F("hnuAccept","simulated neutrino acceptance", 200, 0, 800, 200, -1, 1);
    hnuAccept->GetXaxis()->SetTitle("Electron energy [keV]");
    hnuAccept->GetYaxis()->SetTitle("cos #theta_{#nu}");
    
    TH1D* hNu  = OM.registeredTH1D("hNu","Neutrino spectrum",200,0,800);
    hNu->SetLineColor(3);
    
    TH1F hw("hw","weighting",200,0,1e-29);
   
    size_t npts = 1e9;
    ProgressBar* PB = new ProgressBar(npts, 50);
    for(size_t n=0; n<npts; n++) {
        PB->update(n);
        
        // generate new event
        G.gen_evt_weighted();
        if(G.evt_w <= 0) continue;
       
        double E_e = G.E_2-G.m_2; // electron KE
        // electron-nu cos angle
        double cos_th_enu = G.cos_theta_e_nu();
        
        hSpec->Fill(E_e, G.evt_w);
        hNu->Fill(G.E_1, G.evt_w);
        if(G.K) hw.Fill(G.evt_w);
        
        // vertex position
        RQR.u0[0] = (2*RQR.u0[0]-1)*eCol.r;
        circle_the_square(RQR.u0+1, r_beam);
        
        // momenta, with sign convention reversed
        double p_e[3];
        double p_p[3];
        for(int i=0; i<3; i++) {
            p_e[i] = -G.n_2[i]*G.p_2;
            p_p[i] = -G.p_f[i];
        }
        
        if(p_e[2] > 0) continue; // lose upward-going electrons
        
        pTOF.setInitial(RQR.u0, p_p);
        eTOF.setInitial(RQR.u0, p_e);
        
        double e_passprob = eCol.pass(eTOF);
        if(!e_passprob) continue;
        double passprob = e_passprob * pCol.pass(pTOF);
        double wt = passprob * G.evt_w * G.c_2_wt;
        
        // proton-to-electron time
        double dt = pTOF.calcTOF(pDet_z) - eTOF.calcTOF(eDet_z);
        // wishbone true branch
        bool upper = cos_th_enu < 0;
        
        if(wt) {
            /// uncorrected base asymmetry
            wt *= 1 + a0*G.beta*cos_th_enu;
            
            hWishbone[upper]->Fill(E_e, 1e6*dt, wt);
            haCorn[upper][true]->Fill(E_e, wt);
            
            // non-radiative-corrected component of spectrum:
            // no hard gamma, and evt_w0 (without soft/virtual gamma correction) instead of evt_w
            //if(G.K==0) haCorn[upper][false]->Fill(E_e, G.evt_w0 * G.c_2_wt * passprob);
            
            // [Glu93] parametrized radiative correction
            //haCorn[upper][false]->Fill(E_e, wt*G.Gluck93_radcxn_wt());
            // recoil, weak magnetism correction
            //haCorn[upper][false]->Fill(E_e, wt*G.B59_rwm_cxn_wt());
            
            hPassPos->Fill(RQR.u0[0],RQR.u0[1],wt);
            
            hAccept->Fill(E_e, cos_th_enu, wt);
            hpAccept->Fill(E_e, p_p[2]/G.mag_p_f, wt);
            hnuAccept->Fill(E_e, G.n_1[2], wt);
        }
        
        // re-calculate with grid deflection
        if(true) {
            pTOF.t = pTOF.t_mr;  // time at mirror crossing
            pTOF.calcPos();      // position at mirror crossing
            double Vx = 0, Vy = 0;
            //MD.radial_momentum_deflection(pTOF.xx[0], pTOF.xx[1], Vx, Vy);
            Vx += MD.wires_momentum_deflection(pTOF.xx[0]); // V
            double dpx = Vx / pTOF.v_exit * 2.9979e7; // keV/c
            double dpy = Vy / pTOF.v_exit * 2.9979e7;
            pTOF.kickMomentum(dpx, dpy);
            wt = e_passprob * pCol.pass(pTOF) * G.evt_w * G.c_2_wt;
            wt *= 1 + a0*G.beta*cos_th_enu;
            if(wt) haCorn[upper][false]->Fill(E_e, wt);
        }
    }
    delete PB;
    
    OM.defaultCanvas->SetLeftMargin(0.17);
    OM.defaultCanvas->SetRightMargin(0.04);
    
    hSpec->Draw();
    double acs = 0.8*hSpec->GetMaximum()/haCorn[0][true]->GetMaximum();
    for(int i=0; i<2; i++) {
        haCorn[i][true]->Scale(acs);
        haCorn[i][true]->Draw("Same");
    }
    hNu->Draw("Same");
    OM.printCanvas("Beta_Spectrum");
    
    TH1* hAsym[2];
    int nrebin = 4;
   
    for(int j=0; j<2; j++) {
        hAsym[j] = calcAsym(haCorn[0][j],haCorn[1][j]);
        hAsym[j]->Rebin(nrebin);
        hAsym[j]->Scale(100./nrebin);
        hAsym[j]->SetLineColor(4-2*j);
        hAsym[j]->SetMaximum(0);
        hAsym[j]->GetXaxis()->SetTitle("electron energy [keV]");
        hAsym[j]->GetYaxis()->SetTitle("wishbhone asymmetry [%]");
        hAsym[j]->SetTitle("aCORN simulation");
        hAsym[j]->Draw(j?"HIST Same":"HIST");
        OM.addObject(hAsym[j]);
        printf("Asymmetry %i %.4f\n",j,hAsym[j]->Integral());
    }
    
    TF1 fBeta("fBeta", "[0]*sqrt(x*x+2*511*x)/(511+x)", 100, 400);
    hAsym[true]->Fit(&fBeta, "R+");
    fBeta.SetLineColor(6);
    fBeta.Draw("Same");
    
    OM.printCanvas("Asymmetry");
    
    TF1 lineFit("lineFit", "pol0", 100, 400);
    TH1* hDelta = (TH1*)hAsym[0]->Clone("hDelta");
    hDelta->SetTitle("aCORN near-wires correction");
    hDelta->Add(hAsym[1], -1);
    hDelta->GetYaxis()->SetTitle("false asymmetry [%]");
    hDelta->GetYaxis()->SetTitleOffset(1.4);
    
    //hDelta->Divide(hAsym[1]);
    //hDelta->Scale(100);
    hDelta->SetMinimum(-0.5);
    hDelta->SetMaximum(0.5);
    
    hDelta->Fit(&lineFit, "WWR+");
    hDelta->Draw("HIST");
    lineFit.Draw("Same");
    OM.addObject(hDelta);
    
    // Fred's eLog 135 mirror correction curve
    //TF1 fFredCxn("fFredCxn", "-100*(0.01277 - 4.58e-05*x + 4.13e-08*x*x)", 100, 400);
    //fFredCxn.Draw("Same");
    
    OM.printCanvas("DeltaAsymmetry");
    
    
    gPad->SetCanvasSize(300,300);
    hPassPos->Draw("Col Z");
    OM.printCanvas("Positions");
    
    makeRBpalette();
    
    hAccept->Scale(1.e6/npts);
    hAccept->SetMinimum(-0.5);
    hAccept->SetMaximum(0.5);
    hAccept->Draw("Col");
    OM.printCanvas("Acceptance");
    
    hpAccept->Scale(1.e6/npts);
    hpAccept->SetMinimum(-0.5);
    hpAccept->SetMaximum(0.5);
    hpAccept->Draw("Col");
    OM.printCanvas("pAcceptance");
    
    hnuAccept->Scale(1.e6/npts);
    hnuAccept->SetMinimum(-0.5);
    hnuAccept->SetMaximum(0.5);
    hnuAccept->Draw("Col");
    OM.printCanvas("nuAcceptance");
    
    double wmx = 0.9*hWishbone[0]->GetMaximum();
    hWishbone[1]->Scale(-1.0);
    hWishbone[0]->Add(hWishbone[1]);
    for(int i=0; i<2; i++) {
        hWishbone[i]->SetMinimum(-wmx);
        hWishbone[i]->SetMaximum(wmx);
    }
    hWishbone[0]->Draw("Col");
    //hWishbone[1]->Draw("Col Same");
    OM.printCanvas("Wishbone");
    
    //G.showEffic();
    
    OM.writeROOT();
    
    return 0;
}
