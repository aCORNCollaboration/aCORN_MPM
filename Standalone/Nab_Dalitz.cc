/// \file Nab_Dalitz.cc Electron-proton energy distribution calculation with corrections
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015


// make Nab_Dalitz -j4
// ./Nab_Dalitz

#include "UnpolarizedNeutronDecay.hh"
#include "GraphicsUtils.hh"
#include "GraphUtils.hh"
#include "StringManip.hh"
#include "PathUtils.hh"
#include "OutputManager.hh"
#include "ProgressBar.hh"

#include <TStyle.h>
#include <TGraph.h>
#include <TLatex.h>
#include <gsl/gsl_qrng.h>

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

/// lowest-order Nab distribution given electron kinetic energy [MeV] and proton momentum squared [MeV^2/c^2]
double NabSimple(double KEe, double pp2) {
    if(KEe >= neutronBetaEp/1000.) return 0;
    
    double pe2 = KEe*KEe + 2*(m_e/1000.)*KEe;
    double beta = sqrt(1-m_e*m_e/pow(KEe*1000+m_e,2));
    double Enu = neutronBetaEp/1000. - KEe;
    double f = (pp2-pe2-Enu*Enu) / (2*Enu*sqrt(pe2));
    
    return fabs(f) < 1? (1 + calc_a0()*beta*f) * (KEe+m_e/1000.)*(neutronBetaEp/1000.-KEe) : 0;
}

/// proton p_p^2 as a function of cos theta e_nu, E_e
double pp2(double Ee, double cth) {
    double Enu = delta_mn_mp*1e-3 - Ee;
    double pe2 = Ee*Ee - m_e*m_e*1e-6;
    return 2*cth*sqrt(pe2)*Enu + pe2 + Enu*Enu;
}

/// curve at fixed cos-theta in p_p^2 vs KE_e space
TGraph* dalitz_curve(double cth, int npts = 200) {
    TGraph* g = new TGraph(npts);
    for(int i = 0; i < npts; i++) {
        double Ee = (m_e + neutronBetaEp*double(i)/(npts-1))*1e-3;
        g->SetPoint(i, Ee - m_e*1e-3, pp2(Ee, cth));
    }
    return g;
}

/// make Nab allowed region plot
void draw_Nab_dalitz() {
    TGraph* gm = dalitz_curve(-1);
    TGraph* g0 = dalitz_curve(0);
    TGraph* gp = dalitz_curve(1);
    
    gm->Draw("AC");
    gm->SetMinimum(-0.1);
    gm->SetMaximum(1.5);
    gm->GetXaxis()->SetLimits(-0.05,0.8);
    gm->GetXaxis()->SetTitle("electron kinetic energy [MeV]");
    gm->GetYaxis()->SetTitle("proton momentum p_{p}^{2} [MeV^{2}/c^{2}]");
    gm->SetTitle("");
    
    gm->Draw("AC");
    g0->SetLineStyle(2);
    g0->Draw("C");
    gp->Draw("C"); 
    
    TLatex* Lm = new TLatex(0.42, 0.1, "cos #theta_{e#nu} = -1");
    Lm->Draw();
    
    TLatex* L0 = new TLatex(0.15, 0.5, "cos #theta_{e#nu} = 0");
    L0->Draw();
    
    TLatex* Lp = new TLatex(0.3, 1.37, "cos #theta_{e#nu} = 1");
    Lp->Draw();
}

/// Nab proton distribution at constant E_e
TGraph* Nab_slice(double Ee, int npts = 2) {
    double ppmin = pp2(Ee, -1);
    double ppmax = pp2(Ee, 1);
    TGraph* g = new TGraph(npts + 2);
    g->SetPoint(0, ppmin, 0);
    g->SetPoint(npts+1, ppmax, 0);
    const double gmax = 0.446162;
    for(int i=0; i<npts; i++) {
        double pp = ppmin*1.001 + double(i)/(npts-1)*(ppmax*0.999-ppmin*1.001);
        g->SetPoint(i+1, pp, NabSimple(Ee - m_e*1e-3, pp)/gmax);
    }
    return g;
}

/// make Nab electron slices plot
void draw_Nab_slices() {
    vector<double> KEes = { 0.1, 0.3, 0.5, 0.6 };
    int n=0;
    for(auto KEe: KEes) {
        TGraph* g = Nab_slice(KEe + m_e*1e-3);
        g->SetLineStyle(1+n);
        if(!n) {
            g->Draw("AL");
            g->SetMinimum(0);
            g->SetMaximum(1.1);
            g->GetXaxis()->SetLimits(-0.1, 1.5);
            g->GetXaxis()->SetTitle("proton momentum p_{p}^{2} [MeV^{2}/c^{2}]");
            g->GetXaxis()->SetTitleOffset(1.1);
            g->GetYaxis()->SetTitle("event rate [arb. units]");
            g->SetTitle("");
        } else g->Draw("L");
        n++;
    }
    
    (new TLatex(0.6, 0.95, "100 keV"))->Draw();
    (new TLatex(0.2, 0.82, "300 keV"))->Draw();
    (new TLatex(0.5, 0.69, "500 keV"))->Draw();
    (new TLatex(0.6, 0.51, "600 keV"))->Draw();
   
}


void fill_NabSimple(TH2* h) {
    TAxis* Ax = h->GetXaxis();
    TAxis* Ay = h->GetYaxis();
    Int_t bx,by,bz;
    for(int i=0; i<h->GetNcells(); i++) {
        h->GetBinXYZ(i,bx,by,bz);
        h->SetBinContent(i, NabSimple(Ax->GetBinCenter(bx), Ay->GetBinCenter(by)));
    }
    printf("Simple phase space Max value %g\n", h->GetMaximum());
}

/// Fit a(E) from constant-Ee [MeV] slice
double Nab_slice_fit(TH1* h, double Ee) {
    double ppmin = pp2(Ee, -1);
    double ppmax = pp2(Ee, 1);
    TF1 f("linefit", "pol1", ppmin + 0.05*(ppmax-ppmin), ppmax-0.05*(ppmax-ppmin));
    h->Fit(&f, "QR");
    double x = f.GetParameter(0)/f.GetParameter(1);
    double pe2 = Ee*Ee - m_e*m_e*1e-6;
    double Enu = delta_mn_mp*1e-3 - Ee;
    double a = 2*Ee*Enu / (x + pe2 + Enu*Enu);
    printf("a(%g) = %g\n", Ee, a);
    return a;
}


void Nab_Corrected_Slice(OutputManager& OM, Gluck_beta_MC& G, double KEe, TH1*& h0, TH1*& hRecoil, TH1*& hRad, size_t npts = 1e7) {
    
    h0 = OM.registeredTH1F("h0_"+to_str(KEe),"uncorrected decay distribution", 50, 0, 1.5);
    hRecoil = OM.registeredTH1F("hRecoil_"+to_str(KEe),"recoil-corrected distribution", 50, 0, 1.5);
    hRad = OM.registeredTH1F("hRad_"+to_str(KEe),"radiative-corrected distribution", 50, 0, 1.5);
    vector<TH1*> hs = { h0, hRecoil, hRad };
    for(auto h: hs) {
        h->GetYaxis()->SetTitle("event rate [arb.]");
        h->GetXaxis()->SetTitle("proton momentum^{2} [MeV^{2}/c^{2}]");
    }
    
    G.use_KEe = KEe;
    ProgressBar* PB = new ProgressBar(npts, 50);
    for(size_t n=0; n<npts; n++) {
        PB->update(n);
        
        G.gen_evt_weighted();
        if(G.evt_w <= 0) continue;
        double pp2 = pow(G.mag_p_f/1000, 2);
        
        hRad->Fill(pp2, G.evt_w);
        if(G.K==0) {
            h0->Fill(pp2, G.evt_w0);
            hRecoil->Fill(pp2, G.evt_w0 * G.B59_rwm_cxn_wt());
        }
    }
    delete PB;
    G.use_KEe = -1;
    
    for(auto h: hs) normalize_to_bin_width(h, 1./npts);
    double a0 = Nab_slice_fit(h0, (KEe+m_e)*1e-3);
    double arec = Nab_slice_fit(hRecoil, (KEe+m_e)*1e-3);
    double arad = Nab_slice_fit(hRad, (KEe+m_e)*1e-3);
    printf("da_recoil = %g\nda_rad = %g\n", arec-a0, arad-a0);
    
    hRecoil->Add(h0,-1);
}

int main(int, char**) {
    
    setupSlideStyle();
    gStyle->SetNumberContours(99);
    
    OutputManager OM("Simulated", getEnvSafe("ACORN_SUMMARY")+"/Nab/");
    
    /*
    draw_Nab_dalitz();
    OM.printCanvas("gDalitz");
    draw_Nab_slices();
    OM.printCanvas("gSlices");
    */
    
    GSLQRandom RQR(8);
    Gluck_beta_MC G(&RQR);
    //G.test_calc_P_H();
    G.P_H.q = 0.2;      // more stats for hard gamma component
    
    vector<double> KEes = { 100, 300, 500, 600 };
    vector<TH1*> h0(KEes.size());
    vector<TH1*> hRecoil(KEes.size());
    vector<TH1*> hRad(KEes.size());
    for(size_t n=0; n<KEes.size(); n++) Nab_Corrected_Slice(OM,G, KEes[n], h0[n], hRecoil[n], hRad[n]);
    
    for(size_t n=0; n<KEes.size(); n++) {
        //hRecoil[n]->SetMaximum(1000);
        hRecoil[n]->Draw(n?"HIST SAME":"HIST");
    }
    OM.printCanvas("Nab_recoil");
    
    return EXIT_SUCCESS;
    
    TH2F* hDalitz[2];
    hDalitz[0] = OM.registeredTH2F("hDalitz_0","uncorrected decay distribution", 200, 0, 0.8, 200, 0, 1.5);
    hDalitz[1] = OM.registeredTH2F("hDalitz_1","corrected decay distribution", 200, 0, 0.8, 200, 0, 1.5);
    for(int i=0; i<2; i++) {
        hDalitz[i]->GetXaxis()->SetTitle("Electron energy [MeV]");
        hDalitz[i]->GetYaxis()->SetTitle("proton momentum^{2} [MeV^{2}/c^{2}]");
    }
    
    size_t npts = 1e6;
    ProgressBar* PB = new ProgressBar(npts, 50);
    for(size_t n=0; n<npts; n++) {
        PB->update(n);
        
        G.gen_evt_weighted();
        if(G.evt_w <= 0) continue;
        
        double KEe = (G.E_2 - G.m_2)/1000;
        double pp2 = pow(G.mag_p_f/1000, 2);
        
        // radiative
        //hDalitz[true]->Fill(KEe, pp2, G.evt_w);
        //if(G.K==0) hDalitz[false]->Fill(KEe, pp2, G.evt_w0);
        
        // recoil-order
        if(G.K==0) {
            hDalitz[false]->Fill(KEe, pp2, G.evt_w0);
            hDalitz[true]->Fill(KEe, pp2, G.evt_w0 * G.B59_rwm_cxn_wt());
        }
    }
    delete PB;
    
    fill_NabSimple(hDalitz[0]);
    
    OM.defaultCanvas->SetRightMargin(0.14);
    OM.defaultCanvas->SetLeftMargin(0.13);
    hDalitz[0]->Scale(1./hDalitz[0]->GetMaximum());
    hDalitz[0]->Draw("Col Z");
    OM.printCanvas("hDalitz");
    
    hDalitz[1]->Scale(hDalitz[0]->Integral() / hDalitz[1]->Integral());
    hDalitz[1]->Add(hDalitz[0], -1.0);
    hDalitz[1]->Scale(1000);
    hDalitz[1]->SetMinimum(-4);
    hDalitz[1]->SetMaximum(4);
    makeRBpalette();
    hDalitz[1]->Draw("Col Z");
    OM.printCanvas("hDelta");
}
