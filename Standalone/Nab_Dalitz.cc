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
#include <TGraph2D.h>
#include <TF2.h>
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

TGraph2D* Gluck_r() {
    auto ls = split(loadFileString(getEnvSafe("HOME")+"/Documents/aCORN/20140910_TheoryCorrections_DA100/Gluck_93_Table_II.txt"),"\n");
    vector<double> zs;
    vector< vector<double> > rs;
    vector<double> xs = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    for(auto& l: ls) {
        l = strip(l);
        if(!l.size() || l[0]=='#') continue;
        rs.push_back(sToDoubles(l));
        zs.push_back(rs.back().back());
        rs.back().pop_back();
    }
    
    TGraph2D* g = new TGraph2D(zs.size() * rs[0].size());
    int i = 0;
    for(size_t n=0; n<zs.size(); n++) {
        for(size_t m=0; m<rs[n].size(); m++) {
            g->SetPoint(i++, xs[m], zs[n], rs[n][m]);
        }
    }
    return g;
}

TGraph2D* Gluck_renu() {
    auto ls = split(loadFileString(getEnvSafe("HOME")+"/Documents/aCORN/20140910_TheoryCorrections_DA100/Gluck_93_Table_V.txt"),"\n");
    vector<double> cs;
    vector< vector<double> > rs;
    vector<double> xs = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    for(auto& l: ls) {
        l = strip(l);
        if(!l.size() || l[0]=='#') continue;
        rs.push_back(sToDoubles(l));
        cs.push_back(rs.back().back());
        rs.back().pop_back();
    }
    
    TGraph2D* g = new TGraph2D(cs.size() * rs[0].size());
    int i = 0;
    for(size_t n=0; n<cs.size(); n++) {
        for(size_t m=0; m<rs[n].size(); m++) {
            g->SetPoint(i++, xs[m], cs[n], rs[n][m]);
        }
    }
    
    TF2 fbilinear("fbilinear","[0] + [1]*x + [2]*y + [3]*x*y",0,1,-1,1);
    g->Fit(&fbilinear, "R");
    for(int i=0; i<g->GetN(); i++) {
        double x = g->GetX()[i];
        double c = g->GetY()[i];
        double r = g->GetZ()[i];
        printf("%g\t%g\t%g\t%g\n",x,c,r,100*(r-fbilinear.Eval(x,c)));
    }
    return g;
}

/// Nab proton distribution at constant E_e [MeV]
TGraph* Nab_slice(double Ee, bool with_RWM = false, bool with_RC = false, int npts = 100) {
    static TGraph2D* gGr = Gluck_r();
    double ppmin = pp2(Ee, -1);
    double ppmax = pp2(Ee, 1);
    TGraph* g = new TGraph(npts + 4);
    g->SetPoint(0, -0.1, 0);
    g->SetPoint(1, ppmin, 0);
    g->SetPoint(npts+2, ppmax, 0);
    g->SetPoint(npts+3, 1.5, 0);
    const double gmax = 0.446162;
    double ssum = 0;
    for(int i=0; i<npts; i++) {
        double cos_thn = -0.999 + 1.998*double(i)/(npts-1);
        double pp = ppmin*(1-cos_thn)/2. + ppmax*(cos_thn+1)/2.;
        double s = NabSimple(Ee - m_e*1e-3, pp)/gmax;
        if(with_RWM) s *=  B59_rwm_cxn(Ee*1000, cos_thn);
        if(with_RC) {
            if(cos_thn < -0.975) cos_thn = -0.975;
            if(cos_thn > 0.975) cos_thn = 0.975;
            s *= (1 + 0.01*gGr->Interpolate((Ee*1e3-m_e)/(delta_mn_mp-m_e), 0.5*(cos_thn+1)) + Wilkinson_g_a2pi(Ee*1e3/m_e));
        }
        g->SetPoint(i+2, pp, s);
        ssum += s;
    }
    //printf("ssum = %g\n", ssum);
    return g;
}

/// make Nab electron slices plot
void draw_Nab_slices(bool with_RWM = false, bool with_RC = false) {
    vector<double> KEes = { 0.1, 0.3, 0.5, 0.6 };
    int n=0;
    for(auto KEe: KEes) {
        TGraph* g = Nab_slice(KEe + m_e*1e-3, with_RWM, with_RC);
        g->SetLineStyle(1+n);
        if(!n) {
            g->Draw("AL");
            g->SetMinimum(with_RWM? -1 : 0);
            g->SetMaximum(with_RWM || with_RC? (with_RC? 20:5): 1.1);
            g->GetXaxis()->SetLimits(-0.1, 1.5);
            g->GetXaxis()->SetTitle("proton momentum p_{p}^{2} [MeV^{2}/c^{2}]");
            g->GetXaxis()->SetTitleOffset(1.1);
            g->GetYaxis()->SetTitle(with_RWM || with_RC? "rate difference #times 10^{3}" : "event rate [arb. units]");
            if(with_RC) g->GetYaxis()->SetTitleOffset(1.3);
            g->SetTitle("");
        } else g->Draw("L");
        n++;
    }
    
    if(with_RWM) {
        (new TLatex(0.4, -0.7, "100 keV"))->Draw();
        (new TLatex(-0.05, 3.6, "300 keV"))->Draw();
        (new TLatex(0.42, 4.6, "500 keV"))->Draw();
        (new TLatex(1.05, 2, "600 keV"))->Draw();
    } else if(with_RC) {
        (new TLatex(0.6, 17.0, "100 keV"))->Draw();
        (new TLatex(0.4, 12.5, "300 keV"))->Draw();
        (new TLatex(0.5, 8.7, "500 keV"))->Draw();
        (new TLatex(0.6, 5.0, "600 keV"))->Draw();
    } else {
        (new TLatex(0.6, 0.95, "100 keV"))->Draw();
        (new TLatex(0.2, 0.82, "300 keV"))->Draw();
        (new TLatex(0.5, 0.69, "500 keV"))->Draw();
        (new TLatex(0.6, 0.51, "600 keV"))->Draw();
    }
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
    TF1* f = new TF1(("linefit_"+to_str(Ee)).c_str(), "pol1", ppmin + 0.1*(ppmax-ppmin), ppmax-0.1*(ppmax-ppmin));
    h->Fit(f, "QR");
    double x = f->GetParameter(0)/f->GetParameter(1);
    double pe2 = Ee*Ee - m_e*m_e*1e-6;
    double Enu = delta_mn_mp*1e-3 - Ee;
    double a = 2*Ee*Enu / (x + pe2 + Enu*Enu);
    printf("a(%g) = %g\n", Ee-m_e*1e-3, a);
    return a;
}

/// Fit a(E) from constant-Ee [MeV] slice
double Nab_slice_fit(TGraph* g, double Ee) {
    double ppmin = pp2(Ee, -1);
    double ppmax = pp2(Ee, 1);
    TF1 f("linefit", "pol1", ppmin + 0.1*(ppmax-ppmin), ppmax-0.1*(ppmax-ppmin));
    g->Fit(&f, "QR");
    double x = f.GetParameter(0)/f.GetParameter(1);
    double pe2 = Ee*Ee - m_e*m_e*1e-6;
    double Enu = delta_mn_mp*1e-3 - Ee;
    double a = 2*Ee*Enu / (x + pe2 + Enu*Enu);
    printf("a(%g) = %g\n", Ee-m_e*1e-3, a);
    return a;
}

void setShortPlot(OutputManager& OM) {
    OM.defaultCanvas.SetCanvasSize(200, 100);
    OM.defaultCanvas.SetBottomMargin(0.22);
    OM.defaultCanvas.SetTopMargin(0.03);
    gStyle->SetLabelSize(0.1, "XYZ");
    gStyle->SetTitleSize(0.1,"xyz");
    gStyle->SetTitleOffset(0.7,"y");
    gStyle->SetTitleOffset(0.95,"x");
}

/// Recoil-order or radiative corrections plot
void Nab_cxn_curve(OutputManager& OM, bool radiative = false) {
    int npts = 200;
    TGraph* g = new TGraph(npts);
    for(int n=0; n<npts; n++) {
        double Ee = (m_e + 0.01 + (delta_mn_mp-m_e-0.02)*double(n)/(npts-1))*1e-3;
        TGraph* g0 = Nab_slice(Ee, false);
        TGraph* g1 = Nab_slice(Ee, !radiative, radiative);
        double a0 = Nab_slice_fit(g0,Ee);
        double a1 = Nab_slice_fit(g1,Ee);
        g->SetPoint(n, Ee-m_e*1e-3, 1000*(a1-a0));
    }
    
    TF1 f("linefit", "pol1", 0.05,0.7);
    g->Fit(&f,"N");
    
    setShortPlot(OM);
    
    if(radiative) { 
        g->SetMaximum(1.5);
        g->SetMinimum(-1.5);
        //OM.defaultCanvas.SetLeftMargin(0.19);
    } else g->SetMaximum(0);
    g->Draw("AL");
    g->SetTitle("");
    g->GetXaxis()->SetTitle("electron kinetic energy [MeV]");
    g->GetXaxis()->SetLimits(0,0.8);
    if(radiative) {
        g->GetYaxis()->SetTitle("#delta a^{RC} #times 10^{3}");
        //g->GetYaxis()->SetTitleOffset(0.85);
    }
    else g->GetYaxis()->SetTitle("#delta a^{RWM} #times 10^{3}");
    
    OM.printCanvas(radiative? "Nab_radiative_correction" : "Nab_recoil_correction");
}

/// Gluck simulation-based slice distributions with and without corrections
void Nab_Corrected_Slice(OutputManager& OM, Gluck_beta_MC& G, double KEe, TH1*& h0, TH1*& hRecoil, TH1*& hRad, size_t npts = 1e7) {
    
    h0 = OM.registeredTH1D("h0_"+to_str(KEe),"uncorrected decay distribution", 100, 0, 1.5);
    hRecoil = OM.registeredTH1D("hRecoil_"+to_str(KEe),"recoil-corrected distribution", 100, 0, 1.5);
    hRad = OM.registeredTH1D("hRad_"+to_str(KEe),"radiative-corrected distribution", 100, 0, 1.5);
    vector<TH1*> hs = { h0, hRecoil, hRad };
    for(auto h: hs) {
        h->GetYaxis()->SetTitle("event rate [arb.]");
        h->GetXaxis()->SetTitle("proton momentum p_{p}^{2} [MeV^{2}/c^{2}]");
        h->GetXaxis()->SetTitleOffset(1.1);
    }
    
    G.set_KEe(KEe);
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
    
    for(auto h: hs) normalize_to_bin_width(h, 1./npts);
}

/// Correction cross-check using simulations
void GluckSimCxn(OutputManager& OM) {
    
    vector<double> KEes = { 100, 300, 500, 600 }; // for "Nab_radiative.pdf" individual curves plot
    //vector<double> KEes = { 20, 50, 100, 200, 300, 400, 500, 600, 650, 700, 715, 730, 740, 750, 755, 760 };
    
    vector<TH1*> h0(KEes.size());
    vector<TH1*> hRecoil(KEes.size());
    vector<TH1*> hRad(KEes.size());
    vector<double> da_rad;
    for(size_t n=0; n<KEes.size(); n++) {
        GSLQRandom RQR(8);
        Gluck_beta_MC G(&RQR);
        G.P_H.q = 0.2;      // more stats for hard gamma component
        
        double KEe = KEes[n];
        Nab_Corrected_Slice(OM,G, KEe, h0[n], hRecoil[n], hRad[n], 1e8);
        double a0 = Nab_slice_fit(h0[n], (KEe+m_e)*1e-3);
        double arec = Nab_slice_fit(hRecoil[n], (KEe+m_e)*1e-3);
        double arad = Nab_slice_fit(hRad[n], (KEe+m_e)*1e-3);
        printf("da_recoil = %g\nda_rad = %g\n", arec-a0, arad-a0);
        da_rad.push_back(arad-a0);
    }
    
    if(KEes.size() == 4) { // corrections curves plot
        double nrm = 1.0105*hRecoil[1]->Integral()/hRad[1]->Integral(); // with magic scaling factor to match Gluck '93 table total rate increase
        for(size_t n=0; n<KEes.size(); n++) {
            //hRecoil[n]->SetMaximum(1000);
            //hRecoil[n]->Draw(n?"HIST SAME":"HIST");
            
            hRad[n]->SetLineStyle(1+n);
            hRad[n]->SetLineColor(1);
            hRad[n]->Scale(nrm);
            hRad[n]->Add(h0[n],-1);
            hRad[n]->Scale(1000);
            hRad[n]->SetMinimum(0);
            hRad[n]->SetMaximum(20);
            hRad[n]->SetTitle("");
            hRad[n]->GetYaxis()->SetTitle("rate difference #times 10^{3}");
            hRad[n]->Draw(n?"HIST L SAME":"HIST L");
            (new TLatex(0.6, 17.0, "100 keV"))->Draw();
            (new TLatex(0.4, 12.2, "300 keV"))->Draw();
            (new TLatex(0.5, 8.0, "500 keV"))->Draw();
            (new TLatex(0.6, 4.3, "600 keV"))->Draw();
        }
        OM.printCanvas("Nab_radiative");
    } else { // corrections graph
        setShortPlot(OM);
        TGraph* g = new TGraph(da_rad.size());
        for(size_t n=0; n<KEes.size(); n++) {
            printf("%g\t%g\n", KEes[n], da_rad[n]);
            g->SetPoint(n, KEes[n]/1000., 1000*da_rad[n]);
        }
        
        TF1 f("linefit", "pol1", 0.1,0.7);
        g->Fit(&f,"N");
        
        g->SetMaximum(1.5);
        g->SetMinimum(-1.5);
        g->SetMarkerStyle(7);
        g->Draw("ACP");
        g->SetTitle("");
        g->GetXaxis()->SetTitle("electron kinetic energy [MeV]");
        g->GetXaxis()->SetLimits(0,0.8);
        g->GetYaxis()->SetTitle("#delta a^{RC} #times 10^{3}");
        OM.printCanvas("Nab_MC_radiative_correction");        
    }
}

int main(int, char**) {
    //Gluck_renu(); return 0;
    
    setupSlideStyle();
    gStyle->SetNumberContours(99);
    
    OutputManager OM("Simulated", getEnvSafe("ACORN_SUMMARY")+"/Nab/");
    
    //draw_Nab_dalitz();
    //OM.printCanvas("gDalitz");
    
    //draw_Nab_slices(false,true);
    //OM.printCanvas("Nab_RC_rates");
    
    //Nab_cxn_curve(OM,true);
    
    GluckSimCxn(OM);
    
    return EXIT_SUCCESS;
    
    GSLQRandom RQR(8);
    Gluck_beta_MC G(&RQR);
    //G.test_calc_P_H();
    G.P_H.q = 0.2;      // more stats for hard gamma component
    
    TH2F* hDalitz[2];
    hDalitz[0] = OM.registeredTH2F("hDalitz_0","uncorrected decay distribution", 200, 0, 0.8, 200, 0, 1.5);
    hDalitz[1] = OM.registeredTH2F("hDalitz_1","corrected decay distribution", 200, 0, 0.8, 200, 0, 1.5);
    for(int i=0; i<2; i++) {
        hDalitz[i]->GetXaxis()->SetTitle("electron kinetic energy [MeV]");
        hDalitz[i]->GetYaxis()->SetTitle("proton momentum p_{p}^{2} [MeV^{2}/c^{2}]");
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
    
    OM.defaultCanvas.SetRightMargin(0.14);
    OM.defaultCanvas.SetLeftMargin(0.13);
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
