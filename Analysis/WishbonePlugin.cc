#include "WishbonePlugin.hh"
#include "GraphicsUtils.hh"

WishbonePlugin::WishbonePlugin(RunAccumulator* RA): AnalyzerPlugin(RA,"Wishbone") {
    
    for(uint i=0; i<2; i++) {
        hProtonSignal[i] = registerHist("hProtonSignal_"+itos(i),"Proton detector signal",200,0,20);
        hProtonSignal[i]->GetXaxis()->SetTitle("proton detector ADC channels (#times 10^{3})");
        hProtonSignal[i]->SetLineColor(4-2*i);
    }
    
    hNVeto = registerHist("hNVeto","Veto panels count",9,-0.5,8.5);
    hNVeto->GetXaxis()->SetTitle("Number of veto panels");
    
    unsigned int nEnBins = 500;
    double E0 = 0;
    double E1 = 1000;
    
    TH2F hTemplate("hTemplate","aCORN Wishbone", nEnBins, E0, E1, 500, -0., 10.);
    hTemplate.GetXaxis()->SetTitle("Electron energy [keV]");
    hTemplate.GetYaxis()->SetTitle("Proton TOF [#mus]");
    hWishbone = (TH2F*)registerHist("hWishbone",hTemplate);
    
    TH2F hNETemplate("hNETemplate","PMT trigger counts", 100, 0, 1000, 20, -0.5, 19.5);
    hNETemplate.GetYaxis()->SetTitle("Number of PMTs triggered");
    hNETemplate.GetXaxis()->SetTitle("Electron energy [keV]");
    hNE = (TH2F*)registerHist("hNE",hNETemplate);
}

void WishbonePlugin::fillCoreHists(BaseDataScanner& PDS, double weight) {
    
    if(PDS.E_p_0 <= 0) return;

    hProtonSignal[0]->Fill(PDS.E_p_0/1000., weight);
    
    bool isProton = E_p_lo < PDS.E_p_0 && PDS.E_p_0 < E_p_hi;
    bool isWishboneTime = 3000 < PDS.T_e2p && PDS.T_e2p < 4200;
    
    if(PDS.E_recon > 100 && isWishboneTime) hProtonSignal[1]->Fill(PDS.E_p_0/1000., weight);
    
    if(isProton) {
        hNE->Fill(PDS.E_recon, PDS.nE, weight);
        hNVeto->Fill(PDS.V, weight);
        if(!PDS.V)
            hWishbone->Fill(PDS.E_recon, PDS.T_e2p/1000., weight);
    }
}

void WishbonePlugin::makePlots() {
    
    hWishbone->GetYaxis()->SetRangeUser(2.,5.);
    hWishbone->Draw("Col");
    myA->printCanvas("Wishbone");

    myA->defaultCanvas->SetLogz(true);
    hNE->Draw("Col");
    myA->printCanvas("NPMTs");
    myA->defaultCanvas->SetLogz(false);
    
    myA->defaultCanvas->SetLogy(true);
    
    hNVeto->Draw();
    myA->printCanvas("NVeto");
    
    //hProtonSignal[0]->SetMinimum(hProtonSignal[1]->GetMinimum());
    hProtonSignal[0]->Draw();
    //hProtonSignal[1]->Draw("Same");
    drawVLine(E_p_lo/1000., myA->defaultCanvas, 2);
    drawVLine(E_p_hi/1000., myA->defaultCanvas, 2);
    myA->printCanvas("ProtonSignal");
    
   // myA->defaultCanvas->SetLogy(false);
}
