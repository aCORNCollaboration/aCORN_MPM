#include "aSpectrum.hh"

#include "OutputManager.hh"
#include "GraphUtils.hh"
#include <TStyle.h>
#include <Math/QuasiRandom.h>
#include <TRandom3.h>
#include <TH1F.h>
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

int main(int, char**) {
    
    OutputManager OM("Gluck_MC_Test", "/home/mpmendenhall/Desktop/aCORN_MC/");
    gStyle->SetOptStat("");
    
    TH1F* hBeta_A = new TH1F("hBeta_A", "Analytical beta spectrum", 1000, 0, .800);
    hBeta_A->GetXaxis()->SetTitle("energy [MeV]");
    hBeta_A->GetYaxis()->SetTitle("rate [1/decay/MeV]");
    //TH1F* hBeta_MC = new TH1F("hBeta_MC", "Unweighted MC beta spectrum", 200, 0, .800);
    //hBeta_MC->GetXaxis()->SetTitle("energy [MeV]");
    //hBeta_MC->GetYaxis()->SetTitle("rate [1/decay/MeV]");
    TH1F* hBeta_MC_W = new TH1F("hBeta_MC_W", "Weighted MC beta spectrum", 100, 0, .800);
    hBeta_MC_W->GetXaxis()->SetTitle("energy [MeV]");
    hBeta_MC_W->GetYaxis()->SetTitle("rate [1/decay/MeV]");
    
    // analytical calculation
    double sp = 0;
    BetaSpectrumGenerator BSG(1,1,neutronBetaEp);
    for(int i=1; i<=hBeta_A->GetNbinsX(); i++) {
        double KE = 1000*hBeta_A->GetBinCenter(i);
        double W = (KE+m_e)/m_e;
        //double c_coulomb = WilkinsonF0(1,W); // Gluck MC does not include coulomb correction
        double c_rwm = 1.+Bilenkii59_RWM(W);   // recoil, weak magnetism shape correction
        double c_g = 1.+Wilkinson_g_a2pi(W);   // outer radiative corrections
        double p = plainPhaseSpace(W) * c_g * c_rwm;
        //double p = BSG.decayProb(KE);
        sp += p;
        hBeta_A->SetBinContent(i, p);
    }
    hBeta_A->Rebin(10);
    normalize_to_bin_width(hBeta_A, 1./sp);
    
    RootQRandom RQR;
    Gluck_beta_MC G(&RQR);
    G.test_calc_P_H();
    
    /// weighted simulation mode
    int nToSim = 1e7;
    int nprog = 50;
    double sw = 0;
    printf("Simulating... ");
    for(int i=0; i<nToSim; i++) {
        if(!(i%(nToSim/nprog))) { printf("*"); fflush(stdout); }
        G.gen_evt_weighted();
        if(G.evt_w <= 0) continue;
        G.evt_w *= G.rwm_cxn();
        
        sw += G.evt_w;
        hBeta_MC_W->Fill((G.E_2 - m_e)/1000, G.evt_w);
    }
    printf(" Done.\n");
    normalize_to_bin_width(hBeta_MC_W, 1./sw);
    
    /*
    /// unweighted simulation mode
    printf("Simulating... ");
    sw = 0;
    for(int i=0; i<nToSim; i++) {
        if(!(i%(nToSim/nprog))) { printf("*"); fflush(stdout); }
        G.gen_evt();
        hBeta_MC->Fill((G.E_2 - m_e)/1000);
        sw++;
    }
    printf(" Done.\n");
    normalize_to_bin_width(hBeta_MC, 1./sw);
    */
    
    hBeta_A->SetLineColor(2);
    hBeta_A->Draw("HIST");
    hBeta_MC_W->Draw("HIST Same");
    OM.printCanvas("hBeta");
    
    hBeta_MC_W->Add(hBeta_A, -1.);
    hBeta_MC_W->Draw("HIST");
    OM.printCanvas("hBeta_Resid");
    
    return 0;
}
