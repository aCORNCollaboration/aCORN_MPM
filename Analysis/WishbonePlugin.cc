#include "WishbonePlugin.hh"
#include "GraphicsUtils.hh"
#include "strutils.hh"
#include <TStyle.h>
#include <TColor.h>

int TH2Slicer::pcount = 0;
    
double TH2Slicer::projSlice(double y0, double y1, TH1*& projOut) {
    TAxis* A = h->GetYaxis();
    int b0 = A->FindBin(y0);
    int b1 = A->FindBin(y1);
    
    TH1D* hProj = h->ProjectionX(("_slice_"+itos(pcount++)).c_str(), b0, b1, "e");
    if(projOut) {
        projOut->Add(hProj);
        delete hProj;
    }
    else projOut = hProj;
    
    return A->GetBinUpEdge(b1)-A->GetBinLowEdge(b0);
}

TH2* TH2Slicer::subtractProfile(const TH1* p, double s) const {
    assert(p->GetNbinsX() == h->GetNbinsX());
    TH2* hh = (TH2*)(h->Clone((std::string(h->GetName())+"_bgsub").c_str()));
    for(Int_t nx=0; nx<=h->GetNbinsX()+1; nx++) {
        double dz = p->GetBinContent(nx)*s;
        for(Int_t ny=1; ny<=h->GetNbinsY(); ny++) {
            Int_t b = h->GetBin(nx,ny);
            hh->SetBinContent(b,h->GetBinContent(b)-dz*h->GetYaxis()->GetBinWidth(ny));
        }
    }
    return hh;
}

////////////
////////////
////////////

WishbonePlugin::WishbonePlugin(RunAccumulator* RA): AnalyzerPlugin(RA,"Wishbone") {
    
    for(uint i=0; i<2; i++) {
        hProtonSignal[i] = registerHist("hProtonSignal_"+itos(i),"Proton detector signal",200,0,20);
        hProtonSignal[i]->GetXaxis()->SetTitle("proton detector ADC channels (#times 10^{3})");
        hProtonSignal[i]->SetLineColor(4-2*i);
    }
    
    hNVeto = registerHist("hNVeto","Veto panels count",9,-0.5,8.5);
    hNVeto->GetXaxis()->SetTitle("Number of veto panels");
    
    unsigned int nEnBins = 800;
    double E0 = 0;
    double E1 = 1600;
    
    TH2F hTemplate("hTemplate","aCORN Wishbone", nEnBins, E0, E1, 1001, -0.005, 10.005);
    hTemplate.GetXaxis()->SetTitle("Electron energy [keV]");
    hTemplate.GetYaxis()->SetTitle("Proton TOF [#mus]");
    hTemplate.GetYaxis()->SetTitleOffset(1.4);
    hWishbone = (TH2F*)registerHist("hWishbone",hTemplate);
    hWishbone->GetYaxis()->SetRange();
    
    TH2F hNETemplate("hNETemplate","PMT trigger counts", 100, 0, 1000, 20, -0.5, 19.5);
    hNETemplate.GetYaxis()->SetTitle("Number of PMTs triggered");
    hNETemplate.GetXaxis()->SetTitle("Electron energy [keV]");
    hNE = (TH2F*)registerHist("hNE",hNETemplate);
    
    TH2F hPSTemplate("hPSTemplate","PMT Spectra", 200, 0, 20., NCH_MAX, -0.5, NCH_MAX-0.5);
    hPSTemplate.GetXaxis()->SetTitle("raw signal (#times 10^{3})");
    hPSTemplate.GetYaxis()->SetTitle("detector channel number");
    hPSTemplate.GetYaxis()->SetTitleOffset(1.4);
    hChanSpec = (TH2F*)registerHist("hChanSpec",hPSTemplate);
    
    TH2F hMultTemplate("hMultTemplate","Module Multiplicity", CHAN_PER_MOD+1, -0.5, CHAN_PER_MOD+0.5, CHAN_PER_MOD+1, -0.5, CHAN_PER_MOD+0.5);
    hMultTemplate.GetXaxis()->SetTitle("Module 1 Triggers");
    hMultTemplate.GetYaxis()->SetTitle("Module 2 Triggers");
    hModuleMult = (TH2F*)registerHist("hModuleMult",hMultTemplate);
}

void WishbonePlugin::fillCoreHists(BaseDataScanner& PDS, double weight) {
    
    if(PDS.E_p_0 <= 0) return;
    if(PDS.E_p > 0) hProtonSignal[0]->Fill(PDS.E_p_0/1000., weight);
    
    bool isProton = E_p_lo < PDS.E_p_0 && PDS.E_p_0 < E_p_hi;
    bool isWishboneTime = T_p_lo < PDS.T_e2p && PDS.T_e2p < T_p_hi;
    
    if(PDS.E_p > 0 && PDS.E_recon > 100 && isWishboneTime) hProtonSignal[1]->Fill(PDS.E_p_0/1000., weight);
    
    if(isProton) {
        if(PDS.nE || PDS.nV)
            hModuleMult->Fill(PDS.nFiredMod[0], PDS.nFiredMod[1], weight);
        if(PDS.modDropoutEvt) return;
        for(unsigned int i=0; i<NCH_MAX; i++) {
            if(isWishboneTime && PDS.E_PMT[i])
                hChanSpec->Fill(PDS.E_PMT[i]/1000., i, weight);
        }
        hNE->Fill(PDS.E_recon, PDS.nE, weight);
        hNVeto->Fill(PDS.nV, weight);
        if(!PDS.nV && PDS.T_e2p>0)
            hWishbone->Fill(PDS.E_recon, PDS.T_e2p/1000., weight);
    }
}

void WishbonePlugin::calculateResults() {
    hWishboneEProj[false] = hWishboneEProj[true] = NULL;
    TH2Slicer wbs(hWishbone);
    double tfg = wbs.projSlice(T_p_lo/1000., T_p_hi/1000., hWishboneEProj[true]);
    double tbg = wbs.projSlice(T_p_min/1000., T_p_lo/1000., hWishboneEProj[false]);
    tbg += wbs.projSlice(T_p_hi/1000., T_p_max/1000., hWishboneEProj[false]);
    
    hWishboneBGSub = wbs.subtractProfile(hWishboneEProj[false],1./tbg);
    
    double s0 = 1000./myA->runTimes.total()/hWishboneEProj[true]->GetBinWidth(1);
    hWishboneEProj[true]->Scale(s0);
    hWishboneEProj[false]->Scale(s0*tfg/tbg);
    hWishboneEProj[true]->Add(hWishboneEProj[false],-1.0);
    
    hWishboneEProj[true]->GetYaxis()->SetTitle("rate [mHz/keV]");
    hWishboneEProj[true]->GetYaxis()->SetTitleOffset(1.45);
    hWishboneEProj[true]->SetTitle("aCORN electron spectrum");
    hWishboneEProj[false]->SetTitle("aCORN electron spectrum background");
    hWishboneEProj[false]->SetLineColor(4);
    hWishboneEProj[true]->SetLineColor(2);
    
    hWishboneTProj = hWishbone->ProjectionY("_tproj", 0, -1, "e");
    hWishboneTProj->SetTitle("Wishbone time of flight projection");
    hWishboneTProj->Scale(1./myA->runTimes.total()/hWishboneTProj->GetBinWidth(1));
    hWishboneTProj->GetYaxis()->SetTitle("event rate [Hz/#mus]");
    hWishboneTProj->GetYaxis()->SetTitleOffset(1.45);
    
    hWishboneBGSub->Scale(s0/hWishboneTProj->GetBinWidth(1));
}

void set_plot_style() {
    const Int_t NRGBs = 5;
    const Int_t NCont = 127;
    
    Double_t stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 1.00, 0.75, 1.00 };
    Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.00, 0.80 };
    Double_t blue[NRGBs]  = { 1.00, 0.50, 1.00, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

void WishbonePlugin::makePlots() {
    
    hWishboneBGSub->Rebin2D(4,2);
    hWishboneBGSub->Scale(1./8.);
    hWishboneBGSub->SetMinimum(-5.);
    hWishboneBGSub->SetMaximum(5.);
    set_plot_style();
    hWishboneBGSub->Draw("Col Z");
    myA->printCanvas("Wishbone");
    gStyle->SetPalette(1);
    
    myA->defaultCanvas->SetLogz(true);
    
    hNE->Draw("Col");
    myA->printCanvas("NPMTs");

    hModuleMult->Draw("Col Z");
    myA->printCanvas("ModMult");
    
    hChanSpec->Draw("Col");
    myA->printCanvas("ChannelSpectra");
    
    myA->defaultCanvas->SetLogz(false);
    
    myA->defaultCanvas->SetLogy(true);
    
    TH1* hNVetoR = myA->hToRate(hNVeto,false);
    hNVetoR->GetYaxis()->SetTitle("rate [Hz]");
    hNVetoR->GetYaxis()->SetTitleOffset(1.4);
    hNVetoR->SetMinimum(1e-6);
    hNVetoR->SetMaximum(100);
    hNVetoR->Draw();
    myA->printCanvas("NVeto");
    delete hNVetoR;
    
    //hProtonSignal[0]->SetMinimum(hProtonSignal[1]->GetMinimum());
    TH1* hProtonSignalR = myA->hToRate(hProtonSignal[0],true);
    hProtonSignalR->GetYaxis()->SetTitle("rate [mHz/channel]");
    hProtonSignalR->GetYaxis()->SetTitleOffset(1.4);
    hProtonSignalR->SetMinimum(0.1);
    hProtonSignalR->SetMaximum(100);
    hProtonSignalR->Draw();
    //hProtonSignal[1]->Draw("Same");
    drawVLine(E_p_lo/1000., myA->defaultCanvas, 2);
    drawVLine(E_p_hi/1000., myA->defaultCanvas, 2);
    myA->printCanvas("ProtonSignal");
    delete hProtonSignalR;
    
    myA->defaultCanvas->SetLogy(false);
    int nrebin = 4;
    hWishboneEProj[true]->Rebin(nrebin);
    hWishboneEProj[true]->Scale(1./nrebin);
    hWishboneEProj[false]->Rebin(nrebin);
    hWishboneEProj[false]->Scale(1./nrebin);
    hWishboneEProj[true]->SetMinimum(-0.2);
    hWishboneEProj[true]->SetMaximum(2);
    hWishboneEProj[true]->Draw();
    hWishboneEProj[false]->Draw("Same");
    drawHLine(0., myA->defaultCanvas, 1);
    myA->printCanvas("WishboneEnergy");
    
    hWishboneTProj->SetMinimum(0);
    hWishboneTProj->SetMaximum(5);
    hWishboneTProj->Draw("E0");
    drawVLine(T_p_lo/1000., myA->defaultCanvas, 2);
    drawVLine(T_p_hi/1000., myA->defaultCanvas, 2);
    drawVLine(T_p_min/1000., myA->defaultCanvas, 4);
    drawVLine(T_p_max/1000., myA->defaultCanvas, 4);
    myA->printCanvas("WishboneTime");
}
