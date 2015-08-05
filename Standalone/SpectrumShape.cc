/// \file SpectrumShape.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "aSpectrum.hh"
#include "Collimator.hh"
#include "GraphicsUtils.hh"
#include "StringManip.hh"
#include "PathUtils.hh"
#include "OutputManager.hh"

#include <TStyle.h>
#include <Math/QuasiRandom.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TPad.h>
#include <cassert>

using namespace ROOT::Math;

class RootRandom: public Gluck_MC_Rndm_Src {
public:
    RootRandom() { }
    virtual void next() { R.RndmArray(8,u); }
    TRandom3 R;
};

class RootQRandom: public Gluck_MC_Rndm_Src {
public:
    RootQRandom(): QR_8(11) { }
    virtual void next() { assert(QR_8.Next(u0)); }
    QuasiRandomSobol QR_8;
    double b;
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

int main(int, char**) {
    
    OutputManager OM("Simulated", getEnvSafe("ACORN_SUMMARY")+"/Simulated/");
    
    gStyle->SetOptStat("");
    
    RootQRandom* RQR = new RootQRandom();
    Gluck_beta_MC G(RQR);
    G.test_calc_P_H();
    
    const double B0 = 364;  // field [Gauss]
    const double r_beam = 3.0;
    ElectronTOF eTOF;
    ProtonTOF pTOF;
    CircleCollimator eCol(B0, 2.75, 5.08, 0);
    CircleCollimator pCol(B0, 4.0, 5.08, 0);
    
    G.pt2_max = eCol.pt_max();
    
    ////////////// Wishbone timing estimates
    const double zmin[3] = {0,0,-r_beam};
    const double zmax[3] = {0,0,r_beam};
    const double pemax[3] = {0, 0, sqrt(G.Delta*G.Delta - m_e*m_e)};
    const double pnumax[3] = {0, 0, -(G.Delta-m_e)};
    pTOF.setVertex(zmin, pnumax);
    const double tof_max = pTOF.calcTOF(pDet_z);
    pTOF.setVertex(zmax, pemax);
    const double tof_min = pTOF.calcTOF(pDet_z);
    double tof_mid = 0.5*(tof_max+tof_min);
    printf("Min proton TOF %.2f us, mid %.2f us, max %.2f us\n", tof_min*1e6, tof_mid*1e6, tof_max*1e6);
    tof_mid = 3.4e-6;
    
    int npts = 1e7;
    
    TH1F* hSpec = OM.registeredTH1F("hSpec","Corrected beta spectrum",200,0,800);
    
    TH1F* haCorn[2][2];
    for(int i=0; i<2; i++) {
        for(int j=0; j<2; j++) {
            haCorn[i][j] = OM.registeredTH1F("haCorn_"+to_str(i)+"_"+to_str(j), "asymmetry beta spectrum", 200,0,800);
            haCorn[i][j]->SetLineColor(2+2*i);
        }
    }
    
    TH2F* hPassPos = OM.registeredTH2F("hPassPos","Vertex passed collimator",50,-4,4,50,-4,4);
    
    TH2F* hWishbone = OM.registeredTH2F("hWishbone","simulated wishbone", 400, 0, 800, 400, 2, 5);
    hWishbone->GetXaxis()->SetTitle("Electron energy [keV]");
    hWishbone->GetYaxis()->SetTitle("Proton TOF [#mus]");
    hWishbone->GetYaxis()->SetTitleOffset(1.4);
    
    TH2F* hAccept = OM.registeredTH2F("hAccept","simulated acceptance", 200, 0, 800, 200, -1, 1);
    hAccept->GetXaxis()->SetTitle("Electron energy [keV]");
    hAccept->GetYaxis()->SetTitle("cos #theta_{e#nu}");
    hAccept->GetYaxis()->SetTitleOffset(1.4);
    
    TH2F* hpAccept = OM.registeredTH2F("hpAccept","simulated proton acceptance", 200, 0, 800, 200, -1, 1);
    hpAccept->GetXaxis()->SetTitle("Electron energy [keV]");
    hpAccept->GetYaxis()->SetTitle("cos #theta_{p}");
    hpAccept->GetYaxis()->SetTitleOffset(1.4);
    
    TH2F* hnuAccept = OM.registeredTH2F("hnuAccept","simulated neutrino acceptance", 200, 0, 800, 200, -1, 1);
    hnuAccept->GetXaxis()->SetTitle("Electron energy [keV]");
    hnuAccept->GetYaxis()->SetTitle("cos #theta_{#nu}");
    hnuAccept->GetYaxis()->SetTitleOffset(1.4);
    
    TH1F* hNu  = OM.registeredTH1F("hNu","Neutrino spectrum",200,0,800);
    hNu->SetLineColor(3);
    
    TH1F hw("hw","weighting",200,0,1e-29);
   
    int nprog = 50;
    for(int i=0; i<nprog; i++) printf("="); printf("\n");
    for(int i=0; i<npts; i++) {
        if(!(i%(npts/nprog))) { printf("*"); fflush(stdout); }
        G.gen_evt_weighted();
        if(G.evt_w <= 0) continue;
       
        double E_e = G.E_2-G.m_2; // electron KE
        double cos_th_enu = 0;
        for(int i=0; i<3; i++) cos_th_enu += G.n_1[i]*G.n_2[i]; // electron-nu cos angle
        
        hSpec->Fill(E_e, G.evt_w);
        hNu->Fill(G.E_1, G.evt_w);
        if(G.K) hw.Fill(G.evt_w);
        
        // vertex position
        RQR->u0[0] = (2*RQR->u0[0]-1)*eCol.r;
        circle_the_square(RQR->u0+1, r_beam);
        // momenta, with sign convention reversed
        double p_e[3];
        double p_p[3];
        for(int i=0; i<3; i++) {
            p_e[i] = -G.n_2[i]*G.p_2;
            p_p[i] = -G.p_f[i];
        }
        
        if(p_e[2] > 0) continue; // lose upward-going electrons
        
        double passprob = eCol.pass(RQR->u0, p_e) * pCol.pass(RQR->u0, p_p);
        double wt = passprob * G.evt_w * G.c_2_wt;
        if(!wt) continue;
        
        // proton-to-electron time
        pTOF.setVertex(RQR->u0, p_p);
        eTOF.setVertex(RQR->u0, p_e);
        double dt = pTOF.calcTOF(pDet_z) - eTOF.calcTOF(eDet_z);
        hWishbone->Fill(E_e, 1e6*dt, wt);
        
        bool upper = dt > tof_mid + (tof_min - tof_mid)*(E_e/neutronBetaEp);
        haCorn[upper][true]->Fill(E_e, wt);
        if(G.K==0) haCorn[upper][false]->Fill(E_e, G.evt_w0 * G.c_2_wt * passprob);
        hPassPos->Fill(RQR->u0[0],RQR->u0[1],wt);

        hAccept->Fill(E_e, cos_th_enu, wt);
        hpAccept->Fill(E_e, p_p[2]/G.mag_p_f, wt);
        hnuAccept->Fill(E_e, G.n_1[2], wt);
    }
    printf("\n");
        
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
        hAsym[j]->Scale(1./nrebin);
        hAsym[j]->SetLineColor(4-2*j);
        hAsym[j]->SetMaximum(0);
        hAsym[j]->Draw(j?"HIST Same":"HIST");
        printf("Asymmetry %i %.4f\n",j,hAsym[j]->Integral());
    }
    OM.printCanvas("Asymmetry");
    
    gPad->SetCanvasSize(300,300);
    hPassPos->Draw("Col");
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
    
    double wmx = 0.9*hWishbone->GetMaximum();
    hWishbone->SetMinimum(-wmx);
    hWishbone->SetMaximum(wmx);
    hWishbone->Draw("Col");
    OM.printCanvas("Wishbone");
    
    G.showEffic();
    
    OM.setWriteRoot(true);
    
    return 0;
}
