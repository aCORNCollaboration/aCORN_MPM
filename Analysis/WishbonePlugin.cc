#include "WishbonePlugin.hh"

WishbonePlugin::WishbonePlugin(RunAccumulator* RA): AnalyzerPlugin(RA,"Wishbone") {
    
    hProtonSignal = registerHist("hProtonSignal","Proton detector signal",200,0,20);
    hProtonSignal->GetXaxis()->SetTitle("proton detector ADC channels (#times 10^{3})");
    
    unsigned int nEnBins = 500;
    double E0 = 0;
    double E1 = 1000;
    TH2F hTemplate("hTemplate","aCORN Wishbone", nEnBins, E0, E1, 300, -1., 2.);
    hWishbone = (TH2F*)registerHist("hWishbone",hTemplate);
}

void WishbonePlugin::fillCoreHists(BaseDataScanner& PDS, double weight) {
    hProtonSignal->Fill(PDS.E_p/1000., weight);
    hWishbone->Fill(PDS.E_recon, PDS.T_e2p, weight);
}

void WishbonePlugin::makePlots() {
    hWishbone->Draw();
    myA->printCanvas("Wishbone.pdf");
    
    myA->defaultCanvas->SetLogy(true);
    
    hProtonSignal->Draw();
    myA->printCanvas("ProtonSignal.pdf");
    
    myA->defaultCanvas->SetLogy(false);
}
