#include "aSpectrum.hh"
#include "Collimator.hh"
#include "GraphicsUtils.hh"
#include "strutils.hh"
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
    virtual double selectBranch() { return R.Rndm(); }
    virtual void next_0() { R.RndmArray(5,u); }
    virtual void next_H() { R.RndmArray(8,u); }
    TRandom3 R;
};

class RootQRandom: public Gluck_MC_Rndm_Src {
public:
    RootQRandom(): QR_1(1), QR_5(8), QR_8(11) { }
    virtual double selectBranch() { assert(QR_1.Next(&b)); return b; }
    virtual void next_0() { assert(QR_5.Next(u0)); }
    virtual void next_H() { assert(QR_8.Next(u0)); }
    QuasiRandomSobol QR_1, QR_5, QR_8;
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
    
    RootRandom RR;
    RootQRandom RQR;
    Gluck_beta_MC G(&RQR);
    G.test_calc_P_H();
    
    SimpleCollimator SC;
    ElectronTOF eTOF;
    ProtonTOF pTOF;
    G.pt2_max = SC.pt_max(SC.r_e);
    
    int npts = 1e9;
    
    TH1F* hSpec = OM.registeredTH1F("hSpec","Corrected beta spectrum",200,0,800);
    
    TH1F* haCorn[2][2];
    for(int i=0; i<2; i++) {
        for(int j=0; j<2; j++) {
            haCorn[i][j] = OM.registeredTH1F("haCorn_"+itos(i)+"_"+itos(j), "asymmetry beta spectrum", 200,0,800);
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
   
    for(int i=0; i<npts; i++) {
        if(!(i%(npts/20))) { printf("*"); fflush(stdout); }
        G.gen_evt_weighted();
        if(G.evt_w <= 0) continue;
       
        double E_e = G.E_2-G.m_2; // electron KE
        double cos_th_enu = 0;
        for(int i=0; i<3; i++) cos_th_enu += G.n_1[i]*G.n_2[i]; // electron-nu cos angle
        
        hSpec->Fill(E_e, G.evt_w);
        
        hNu->Fill(G.E_1, G.evt_w);
        
        if(G.K) hw.Fill(G.evt_w);
        
        RQR.u0[0] = 4*RQR.u0[0]-2;
        circle_the_square(RQR.u0+1, 2.0);
        
        for(int i=0; i<3; i++) {
            SC.x[i] = RQR.u0[i];
            SC.p_e[i] = -G.n_2[i]*G.p_2;
            SC.p_p[i] = -G.p_f[i];
        }
        
        if(SC.pass()) {
            haCorn[G.n_1[2] > 0][true]->Fill(E_e, 1000*G.evt_w*G.c_2_wt);
            if(G.K==0) haCorn[G.n_1[2] > 0][false]->Fill(E_e, 1000*G.evt_w0*G.c_2_wt);
            hPassPos->Fill(SC.x[0],SC.x[1],G.evt_w*G.c_2_wt);

            hAccept->Fill(E_e, cos_th_enu, G.evt_w*G.c_2_wt);
            hpAccept->Fill(E_e, SC.p_p[2]/G.mag_p_f, G.evt_w*G.c_2_wt);
            hnuAccept->Fill(E_e, G.n_1[2], G.evt_w*G.c_2_wt);
            
            // proton-to-electron time
            double dt = pTOF.calcTOF(SC.x, SC.p_p) - eTOF.calcTOF(SC.x, SC.p_e);
            hWishbone->Fill(E_e, 1e6*dt, G.evt_w*G.c_2_wt);
        }
    }
    printf("\n");
        
    hSpec->Draw();
    for(int i=0; i<2; i++) haCorn[i][true]->Draw("Same");
    hNu->Draw("Same");
    OM.printCanvas("Beta_Spectrum");
    
    TH1* hAsym[2];
    int nrebin = 4;
    for(int j=0; j<2; j++) {
        hAsym[j] = calcAsym(haCorn[0][j],haCorn[1][j]);
        hAsym[j]->Rebin(nrebin);
        hAsym[j]->Scale(1./nrebin);
        hAsym[j]->SetLineColor(4-2*j);
        hAsym[j]->Draw(j?"Same":"");
        printf("Asymmetry %i %.4f\n",j,hAsym[j]->Integral());
    }
    OM.printCanvas("Asymmetry");
    
    gPad->SetCanvasSize(300,300);
    hPassPos->Draw("Col");
    OM.printCanvas("Positions");
    
    makeRBpalette();
    
    hAccept->Scale(1.e6/npts);
    hAccept->SetMinimum(-0.3);
    hAccept->SetMaximum(0.3);
    hAccept->Draw("Col");
    OM.printCanvas("Acceptance");
    
    hpAccept->Scale(1.e6/npts);
    hpAccept->SetMinimum(-0.3);
    hpAccept->SetMaximum(0.3);
    hpAccept->Draw("Col");
    OM.printCanvas("pAcceptance");
    
    hnuAccept->Scale(1.e6/npts);
    hnuAccept->SetMinimum(-0.3);
    hnuAccept->SetMaximum(0.3);
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
