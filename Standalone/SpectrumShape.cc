#include "aSpectrum.hh"
#include <Math/QuasiRandom.h>
#include <TRandom3.h>
#include <TH1F.h>
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
    RootQRandom(): QR_1(1), QR_5(5), QR_8(8) { }
    virtual double selectBranch() { assert(QR_1.Next(&b)); return b; }
    virtual void next_0() { assert(QR_5.Next(u)); }
    virtual void next_H() { assert(QR_8.Next(u)); }
    QuasiRandomSobol QR_1, QR_5, QR_8;
    double b;
};

int main(int, char**) {
    Gluck_beta_MC G(new RootQRandom());
    G.test_calc_P_H();
    
    int npts = 1e6;
    
    TH1F hSpec("hSpec","Corrected beta spectrum",200,0,800);
    TH1F hSpec2("hSpec2","other beta spectrum",200,0,800);
    hSpec2.SetLineColor(2);
    TH1F hNu("hNu","Neutrino spectrum",200,0,800);
    hNu.SetLineColor(3);
    
    TH1F hw("hw","weighting",200,0,1e-29);
   
    
    for(int i=0; i<npts; i++) {
        if(!(i%(npts/20))) { printf("*"); fflush(stdout); }
        G.gen_evt_weighted();
        hSpec.Fill(G.E_2-G.m_2, G.evt_w);
        hSpec2.Fill(G.E_2-G.m_2, G.evt_w0);
        hNu.Fill(G.E_1, G.evt_w);
        if(G.K) hw.Fill(G.evt_w);
    }
    printf("\n");
    /*
    G.SetRandom(new RootRandom());
    
    for(int i=0; i<npts; i++) {
        if(!(i%(npts/20))) { printf("*"); fflush(stdout); }
        G.gen_evt();
        hSpec2.Fill(G.E_2-G.m_2);
    }
    printf("\n");
    */
    
    hSpec.Draw();
    hSpec2.Draw("Same");
    hNu.Draw("Same");
    gPad->Print("Gluck_Beta_Spectrum.pdf");
    
    gPad->SetLogy(true);
    hw.Draw();
    gPad->Print("Gluck_w_H.pdf");
    
    G.showEffic();
    
    return 0;
}
