#include "WishbonePlugin.hh"

WishbonePlugin::WishbonePlugin(RunAccumulator* RA): AnalyzerPlugin(RA,"Wishbone") {
    
    hProtonSignal = registerHist("hProtonSignal","Proton detector signal",200,0,20);
    hProtonSignal->GetXaxis()->SetTitle("proton detector ADC channels (#times 10^{3})");
    
    unsigned int nEnBins = 500;
    double E0 = 0;
    double E1 = 1000;
    TH2F hTemplate("hTemplate","aCORN Wishbone", nEnBins, E0, E1, 200, -0., 8.);
    hTemplate.GetXaxis()->SetTitle("Electron energy [keV]");
    hTemplate.GetYaxis()->SetTitle("Proton TOF [#mus]");
    hWishbone = (TH2F*)registerHist("hWishbone",hTemplate);
}

void WishbonePlugin::fillCoreHists(BaseDataScanner& PDS, double weight) {
    hProtonSignal->Fill(PDS.E_p/1000., weight);
    if(500 < PDS.E_p && PDS.E_p < 2000)
        hWishbone->Fill(PDS.E_recon, PDS.T_e2p/1000., weight);
}

void WishbonePlugin::makePlots() {
    
    hWishbone->Draw("Col");
    myA->printCanvas("Wishbone");
    
    myA->defaultCanvas->SetLogy(true);
    
    hProtonSignal->Draw();
    myA->printCanvas("ProtonSignal");
    
    myA->defaultCanvas->SetLogy(false);
}
