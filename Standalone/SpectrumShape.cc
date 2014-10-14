#include "aSpectrum.hh"
#include "Collimator.hh"
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
    RootRandom RR;
    RootQRandom RQR;
    Gluck_beta_MC G(&RQR);
    G.test_calc_P_H();
    
    SimpleCollimator SC;
    G.pt2_max = SC.pt_max(SC.r_e);
    
    int npts = 1e8;
    
    TH1F hSpec("hSpec","Corrected beta spectrum",200,0,800);
    
    TH1F haCorn[2] = { TH1F("haCorn1","other beta spectrum",200,0,800), TH1F("haCorn2","other beta spectrum",200,0,800) };
    for(int i=0; i<2; i++) haCorn[i].SetLineColor(2+2*i);
    
    TH2F hPassPos("hPassPos","Vertex passed collimator",50,-4,4,50,-4,4);
    
    TH1F hNu("hNu","Neutrino spectrum",200,0,800);
    hNu.SetLineColor(3);
    
    TH1F hw("hw","weighting",200,0,1e-29);
   
    
    for(int i=0; i<npts; i++) {
        if(!(i%(npts/20))) { printf("*"); fflush(stdout); }
        G.gen_evt_weighted();
        if(G.evt_w <= 0) continue;
        
        hSpec.Fill(G.E_2-G.m_2, G.evt_w);
        //hSpec2.Fill(G.E_2-G.m_2, G.evt_w0);
        
        hNu.Fill(G.E_1, G.evt_w);
        
        if(G.K) hw.Fill(G.evt_w);
        
        for(int i=0; i<3; i++) {
            SC.x[i] = 4*RQR.u0[i]-2;
            SC.p_e[i] = G.n_2[i]*G.p_2;
            SC.p_p[i] = G.p_f[i];
        }
        if(SC.pass()) {
            haCorn[G.n_1[2] > 0].Fill(G.E_2-G.m_2, 1000*G.evt_w*G.c_2_wt);
            hPassPos.Fill(SC.x[0],SC.x[1],G.evt_w*G.c_2_wt);
        }
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
    for(int i=0; i<2; i++) haCorn[i].Draw("Same");
    hNu.Draw("Same");
    gPad->Print("Gluck_Beta_Spectrum.pdf");
    
    TH1* hAsym = calcAsym(&haCorn[0],&haCorn[1]);
    hAsym->Draw();
    gPad->Print("Gluck_Asymmetry.pdf");
    
    gPad->SetCanvasSize(300,300);
    hPassPos.Draw("Col");
    gPad->Print("Gluck_Positions.pdf");
    
    //gPad->SetLogy(true);
    //hw.Draw();
    //gPad->Print("Gluck_w_H.pdf");
    
    G.showEffic();
    
    return 0;
}
