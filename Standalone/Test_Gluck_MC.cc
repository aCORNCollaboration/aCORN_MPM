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

/*
class RootQRandom: public Gluck_MC_Rndm_Src {
public:
    RootQRandom(): QR_8(11) { }
    virtual void next() { assert(QR_8.Next(u0)); }
    QuasiRandomSobol QR_8;
    double b;
};
*/

class RootQRandom: public Gluck_MC_Rndm_Src {
public:
    RootQRandom(): QR_8(8) { }
    virtual void next() { assert(QR_8.Next(u)); }
    QuasiRandomSobol QR_8;
    double b;
};

double dot3x(const double* a, const double* b) {
    double s = 0;
    for(int i=0; i<3; i++) s += a[i]*b[i];
    return s;
}

int main(int, char**) {
    
    OutputManager OM("Gluck_MC_Test", "/home/mpmendenhall/Desktop/aCORN_MC/");
    gStyle->SetOptStat("");
    
    TH1F* hBeta_A[2];
    TH1F* hBeta_MC[2];
    for(int i=0; i<2; i++) {
        hBeta_A[i] = new TH1F(i?"hBeta_A1":"hBeta_A0", "Analytical beta spectrum", 1000, 0, .800);
        hBeta_A[i]->GetXaxis()->SetTitle("energy [MeV]");
        hBeta_A[i]->GetYaxis()->SetTitle("rate [1/decay/MeV]");
        
        hBeta_MC[i] = new TH1F(i?"hBeta_MC1":"hBeta_MC0", "MC beta spectrum", 50, 0, .800);
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
        hBeta_A[i]->Rebin(20);
        normalize_to_bin_width(hBeta_A[i], 1./sp);
    }
    
    RootQRandom RQR;
    Gluck_beta_MC G(&RQR);
    //G.test_calc_P_H();
    
    /// weighted simulation mode
    int nToSim = 1e7;
    int nprog = 50;
    double sw = 0;
    printf("Simulating... ");
    for(int i=0; i<nToSim; i++) {
        if(!(i%(nToSim/nprog))) { printf("*"); fflush(stdout); }
        G.gen_evt_weighted();
        if(G.evt_w <= 0) continue;
        double cos_thn = dot3x(G.n_1, G.n_2);
        G.evt_w *= G.coulomb_cxn()*G.rwm_cxn();
        
        sw += G.evt_w;
        hBeta_MC[cos_thn > 0]->Fill((G.E_2 - m_e)/1000, G.evt_w);
    }
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
    
    return 0;
}
