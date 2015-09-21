/// \file WishbonePlugin.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "WishbonePlugin.hh"
#include "WishboneFit.hh"
#include "PathUtils.hh"
#include "GraphicsUtils.hh"
#include "GraphUtils.hh"
#include "StringManip.hh"
#include <TStyle.h>
#include <TColor.h>

int TH2Slicer::pcount = 0;
    
double TH2Slicer::projSlice(double y0, double y1, TH1*& projOut) {
    TAxis* A = h->GetYaxis();
    int b0 = A->FindBin(y0);
    int b1 = A->FindBin(y1);
    
    TH1D* hProj = h->ProjectionX(("_slice_"+to_str(pcount++)).c_str(), b0, b1, "e");
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

WishbonePlugin::WishbonePlugin(RunAccumulator* RA, const string& nm, const string& inflname):
RunAccumulatorPlugin(RA, nm, inflname),
hProtonSignal(this), hNVeto(this), hVetoSum(this), hNE(this), hPMTs(this),
hChanSpec(this), hModuleMult(this), hPos(this), hPosSigma(this), hEnergyRadius(this) {
    
    dataMode = BAD;
    string dsrc = split(strip(getEnvSafe("ACORN_REDUCED_ROOT"),"/"),"/").back();
    if(dsrc == "ROOT_NG6") { dataMode = NG6; config_NG6_cuts(); }
    if(dsrc == "ROOT_NGC") { dataMode = NGC; config_NGC_cuts(); }
    assert(dataMode != BAD);
    
    TH1F hProtonTemplate("hProtonSignal", "Proton detector signal", 200,0,20);
    hProtonTemplate.GetXaxis()->SetTitle("proton detector ADC channels (#times 10^{3})");
    initRegions(hProtonSignal);
    hProtonSignal.setTemplate(hProtonTemplate);
  
    TH1F hNVetoTemplate("hNVeto","Veto panels count",9,-0.5,8.5);
    hNVetoTemplate.GetXaxis()->SetTitle("Number of veto panels");
    initRegions(hNVeto);
    hNVeto.setTemplate(hNVetoTemplate);
    
    TH1F hVetoSumTemplate("hVetoSum","Veto PMTs signal sum",100,0,20);
    hVetoSumTemplate.GetXaxis()->SetTitle("raw signal sum (#times 10^{3})");
    initRegions(hVetoSum);
    hVetoSum.setTemplate(hVetoSumTemplate);
    
    unsigned int nEnBins = 400;
    double E0 = 0;
    double E1 = 1600;
    double timeBin = 0.02; // time bin width, microseconds
    
    // note 5ns shift in time axis bins, to avoid bin boundaries at 10ns+/-error unit boundaries.
    TH2F hTemplate("hTemplate","aCORN Wishbone", nEnBins, E0, E1, 10/timeBin + 1, - 0.005, 10. + timeBin - 0.005);
    hTemplate.GetXaxis()->SetTitle("Electron energy [keV]");
    hTemplate.GetYaxis()->SetTitle("Proton TOF [#mus]");
    hWishbone = (TH2F*)registerSavedHist("hWishbone",hTemplate);
    hWishbone->SetTitle("aCORN NG-C Wishbone");
    hWishbone->GetYaxis()->SetRange();
    
    TH2F hNETemplate("hNE","PMT trigger counts", 100, 0, 1000, 20, -0.5, 19.5);
    hNETemplate.GetYaxis()->SetTitle("Number of PMTs triggered");
    hNETemplate.GetXaxis()->SetTitle("Electron energy [keV]");
    initRegions(hNE);
    hNE.setTemplate(hNETemplate);
    
    TH2F hPSTemplate("hChanSpec","PMT Spectra", 200, 0, (1<<15)/1000., NCH_MAX, -0.5, NCH_MAX-0.5);
    hPSTemplate.GetXaxis()->SetTitle("raw signal (#times 10^{3})");
    hPSTemplate.GetYaxis()->SetTitle("detector channel number");
    initRegions(hChanSpec);
    hChanSpec.setTemplate(hPSTemplate);
    
    TH2F hMultTemplate("hModuleMult","Module Multiplicity", CHAN_PER_MOD+1, -0.5, CHAN_PER_MOD+0.5, CHAN_PER_MOD+1, -0.5, CHAN_PER_MOD+0.5);
    hMultTemplate.GetXaxis()->SetTitle("Module 1 Triggers");
    hMultTemplate.GetYaxis()->SetTitle("Module 2 Triggers");
    initRegions(hModuleMult);
    hModuleMult.setTemplate(hMultTemplate);
    
    TH2F hPosTemplate("hPos","Beta hit positions", 100, -3, 3, 100, -3, 3);
    hPosTemplate.GetXaxis()->SetTitle("x [PMT spacings]");
    hPosTemplate.GetYaxis()->SetTitle("y [PMT spacings]");
    initRegions(hPos);
    hPos.setTemplate(hPosTemplate);
    
    TH2F hPosSigmaTemplate("hPosSigma","Beta hit position spread", 100, 0, 2, 100, 0, 2);
    hPosSigmaTemplate.GetXaxis()->SetTitle("#sigma_{x} [PMT spacings]");
    hPosSigmaTemplate.GetYaxis()->SetTitle("#sigma_{y} [PMT spacings]");
    initRegions(hPosSigma);
    hPosSigma.setTemplate(hPosSigmaTemplate);
    
    
    TH2F hEnergyRadiusTemplate("hEnergyRadius","Beta hit positions", 100, 0, 4, 100, 0, 800);
    hEnergyRadiusTemplate.GetXaxis()->SetTitle("radius^{2} [(PMT spacings)^{2}]");
    hEnergyRadiusTemplate.GetYaxis()->SetTitle("Energy [keV]");
    initRegions(hEnergyRadius);
    hEnergyRadius.setTemplate(hEnergyRadiusTemplate);
}

void WishbonePlugin::config_NGC_cuts() {
    E_p_lo = 550;
    E_p_hi = 1500;
    T_p_min = 750;
    T_p_lo = 3000;
    T_p_hi = 4500;
    T_p_max = 9500;
}

void WishbonePlugin::config_NG6_cuts() {
    E_p_lo = 650;
    E_p_hi = 2400;
    T_p_min = 750;
    T_p_lo = 2750;
    T_p_hi = 4500;
    T_p_max = 9500;
}
   
void WishbonePlugin::initRegions(FGBGRegionsHist& h) {
    h.addRegion(T_p_min, T_p_lo, false);
    h.addRegion(T_p_lo, T_p_hi, true);
    h.addRegion(T_p_hi, T_p_max, false);
}

void WishbonePlugin::fillCoreHists(BaseDataScanner& PDS, double weight) {
    
    for(unsigned int i=0; i<NCH_MAX; i++) {
        if(PDS.E_PMT[i]) {
            assert(PDS.E_PMT[i] >= 0);
            hChanSpec.fill(3500, PDS.E_PMT[i]/1000., i, weight);
        }
    }
    
    if(PDS.E_p_0 <= 0) return;
    
    hProtonSignal.fill(PDS.T_e2p, PDS.E_p_0/1000., weight);
    
    bool isProton = E_p_lo < PDS.E_p_0 && PDS.E_p_0 < E_p_hi;
    
    if(isProton) {
        if(PDS.nE || PDS.nV)
            hModuleMult.fill(PDS.T_e2p, PDS.nFiredMod[0], PDS.nFiredMod[1], weight);
        if(PDS.modDropoutEvt) return;
        
        double vetosum = 0;
        for(unsigned int i=0; i<NCH_MAX; i++) {
            if(PDS.E_PMT[i]) {
                //hChanSpec.fill(PDS.T_e2p, PDS.E_PMT[i]/1000., i, weight);
                if(i<N_V_PMT) vetosum += PDS.E_PMT[i];
            }
        }
        if(vetosum) hVetoSum.fill(PDS.T_e2p, vetosum/1000., weight);
            
        hNE.fill(PDS.T_e2p, PDS.E_recon, PDS.nE, weight);
        hNVeto.fill(PDS.T_e2p, PDS.nV, weight);
        if(!PDS.nV && PDS.T_e2p>0)
            hWishbone->Fill(PDS.E_recon, PDS.T_e2p/1000., weight);
        
        hPos.fill(PDS.T_e2p, PDS.Pos.px[0], PDS.Pos.px[1], weight);
        hPosSigma.fill(PDS.T_e2p, PDS.Pos.sx[0], PDS.Pos.sx[1], weight);
        hEnergyRadius.fill(PDS.T_e2p, PDS.Pos.r2(), PDS.E_recon, weight);
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
    
    hWishboneEProj[true]->GetYaxis()->SetTitle("rate [Hz/MeV]");
    hWishboneEProj[true]->SetTitle("aCORN electron spectrum");
    hWishboneEProj[false]->SetTitle("aCORN electron spectrum background");
    hWishboneEProj[false]->SetLineColor(4);
    hWishboneEProj[true]->SetLineColor(2);
    
    hWishboneTProj = hWishbone->ProjectionY("_tproj", 0, -1, "e");
    hWishboneTProj->SetTitle("Wishbone time of flight projection");
    normalize_to_bin_width(hWishboneTProj, 1./myA->runTimes.total());
    hWishboneTProj->GetYaxis()->SetTitle("event rate [Hz/#mus]");
    
    hWishboneBGSub->Scale(s0/hWishboneTProj->GetBinWidth(1));
    hWishboneBGSub->GetZaxis()->SetTitle("rate [Hz/MeV/#mus]");
    
    int eb0 = hWishboneBGSub->GetXaxis()->FindBin(100);
    int eb1 = hWishboneBGSub->GetXaxis()->FindBin(300);
    hWishboneFiducialTProj = hWishboneBGSub->ProjectionY("_tproj", eb0, eb1, "e");
    hWishboneFiducialTProj->Scale(hWishboneBGSub->GetXaxis()->GetBinWidth(1)/1000.);
    hWishboneFiducialTProj->GetYaxis()->SetTitle("event rate [Hz/#mus]");

    hWBRate = (TH2F*)hWishbone->Clone("hWBRate");
    normalize_to_bin_area(hWBRate, 1000./myA->runTimes.total());
    hWBRate->GetZaxis()->SetTitle("rate [Hz/MeV/#mus]");
    
    //-------------------
    
    hPos.makeRates(2);
    hPosSigma.makeRates(2);
    hEnergyRadius.makeRates(2);
    hChanSpec.makeRates(2);
    hNE.makeRates(0);
    hModuleMult.makeRates(0);
    hNVeto.makeRates(0);
    hNVeto.hRates[true]->GetYaxis()->SetTitle("rate [Hz]");
    
    hProtonSignal.makeRates(1);
    hProtonSignal.hRates[false]->GetYaxis()->SetTitle("rate [Hz/kchannel]");
   
    hVetoSum.makeRates(1,1000.);
    hVetoSum.hRates[false]->GetYaxis()->SetTitle("rate [#muHz/channel]");
    
    int nrebin = 2;
    for(int i=0; i<2; i++) {
        hWishboneEProj[i]->Rebin(nrebin);
        hWishboneEProj[i]->Scale(1./nrebin);
    }
    hWishboneEProj[false]->GetYaxis()->SetTitle("Background rate [Hz/MeV]");
}

void WishbonePlugin::makeAnaResults() {
    AnaResult baseResult = myA->makeBaseResult();
    
    baseResult.value = integralAndError(hProtonSignal.hRates[false], E_p_lo/1000., E_p_hi/1000., baseResult.err, "width");
    AcornDB::ADB().uploadAnaResult("proton_bg", "Proton peak background rate [Hz]", baseResult);
    baseResult.value = integralAndError(hProtonSignal.hRates[true], E_p_lo/1000., E_p_hi/1000., baseResult.err, "width");
    AcornDB::ADB().uploadAnaResult("proton_fg", "Proton peak foreground rate [Hz]", baseResult);
    
    baseResult.value = integralAndError(hWishboneEProj[false], 1., 1000., baseResult.err, "width")/1000;
    baseResult.err /= 1000;
    AcornDB::ADB().uploadAnaResult("wb_bg", "Wishbone background rate [Hz]", baseResult);
    baseResult.value = integralAndError(hWishboneEProj[true], 1., 1000., baseResult.err, "width")/1000;
    baseResult.err /= 1000;
    AcornDB::ADB().uploadAnaResult("wb_fg", "Wishbone foreground rate [Hz]", baseResult);
    
    baseResult.value = integralAndError(hWishboneFiducialTProj, 3.0, 3.6, baseResult.err, "width");
    AcornDB::ADB().uploadAnaResult("wb_fast_fiducial", "fast protons in energy fiducial [Hz]", baseResult);
    baseResult.value = integralAndError(hWishboneFiducialTProj, 3.7, 4.5, baseResult.err, "width");
    AcornDB::ADB().uploadAnaResult("wb_slow_fiducial", "slow protons in energy fiducial [Hz]", baseResult);
    
    ManualWishboneSeparator MWS("WBFit", myA);
    MWS.setWishbone(hWishboneBGSub);
    MWS.extractAsymmetry();
    baseResult.value = MWS.fitAsym(50,400,baseResult.err);
    AcornDB::ADB().uploadAnaResult("obs_asym", "observed uninterpreted asymmetry", baseResult);
}

void WishbonePlugin::makePlots() {
    
    defaultCanvas->SetRightMargin(0.20);

    bool isCombined = myA->runTimes.nTags() > 1000;
    if(true || !isCombined) {
        hWishboneBGSub->Rebin2D(2,2);
        hWishboneBGSub->Scale(1./4.);
    }
    hWishboneBGSub->SetMinimum(-10.);
    hWishboneBGSub->SetMaximum(10.);
    hWishboneBGSub->GetXaxis()->SetRangeUser(0,1000);
    hWishboneBGSub->GetYaxis()->SetRangeUser(2,5);
    makeRBpalette();
    hWishboneBGSub->Draw("Col Z");
    drawHLine(T_p_lo/1000., defaultCanvas, 2);
    drawHLine(T_p_hi/1000., defaultCanvas, 2);
    drawHLine(T_p_min/1000., defaultCanvas, 4);
    drawHLine(T_p_max/1000., defaultCanvas, 4);
    printCanvas("Wishbone");
    
    defaultCanvas->SetRightMargin(0.15);
    
    makeGrayscalepalette(false);
    hWBRate->SetMaximum(10);
    hWBRate->SetContour(255);
    hWBRate->GetXaxis()->SetRangeUser(0,1000);
    hWBRate->GetYaxis()->SetRangeUser(2,5);
    hWBRate->Draw("Col Z");
    printCanvas("Wishbone_NoSub");
    
    gStyle->SetPalette(1);
      
    hPos.hRates[true]->SetMinimum(0);
    hPos.hRates[true]->SetMaximum(1);
    hPos.hRates[true]->Draw("Col Z");
    Positioner Pos;
    Pos.drawPMTs(1,true);
    printCanvas("Positions");
    
    hPosSigma.hRates[true]->SetMinimum(0);
    hPosSigma.hRates[true]->SetMaximum(10);
    hPosSigma.hRates[true]->Draw("Col Z");
    printCanvas("PosSigma");
    
    hEnergyRadius.hRates[true]->SetMinimum(0);
    hEnergyRadius.hRates[true]->SetMaximum(0.03);
    hEnergyRadius.hRates[true]->Draw("Col Z");
    printCanvas("EnergyRadius");
    
    defaultCanvas->SetLogz(true);
    
    //hChanSpec.hRates[true]->SetMaximum(0.4);
    hChanSpec.hRates[true]->GetXaxis()->SetRangeUser(0,16);
    hChanSpec.hRates[true]->Draw("Col Z");
    printCanvas("ChannelSpectra");    
    hChanSpec.hRates[false]->GetXaxis()->SetRangeUser(0,16);
    hChanSpec.hRates[false]->Draw("Col Z");
    printCanvas("ChannelSpectraBG"); 
    
    hNE.hRates[true]->Draw("Col Z");
    printCanvas("NPMTs");

    hModuleMult.hRates[true]->Draw("Col Z");
    printCanvas("ModMult");
    
    defaultCanvas->SetLogz(false);
    
    defaultCanvas->SetRightMargin(0.04);
    
    defaultCanvas->SetLogy(true);

    hNVeto.hRates[true]->SetMinimum(1e-7);
    hNVeto.hRates[true]->SetMaximum(10);
    hNVeto.hRates[true]->Draw();
    hNVeto.hRates[false]->Draw("Same");
    printCanvas("NVeto");
    
    hProtonSignal.hRates[false]->SetMinimum(isCombined?1e-4:1e-3);
    hProtonSignal.hRates[false]->SetMaximum(50);
    hProtonSignal.hRates[false]->Draw();
    hProtonSignal.hRates[true]->Draw("Same");
    drawVLine(E_p_lo/1000., defaultCanvas, 2);
    drawVLine(E_p_hi/1000., defaultCanvas, 2);
    printCanvas("ProtonSignal");
  
    defaultCanvas->SetLogy(false);
    
    hVetoSum.hRates[false]->SetMinimum(-20);
    hVetoSum.hRates[false]->SetMaximum(200);
    hVetoSum.hRates[false]->Draw();
    drawHLine(0, defaultCanvas, 1);
    hVetoSum.hRates[true]->Draw("Same");
    printCanvas("VetoSum");
    
    hWishboneEProj[true]->SetMinimum(-0.2);
    hWishboneEProj[true]->SetMaximum(dataMode == NG6? 3 : 15);
    hWishboneEProj[true]->GetXaxis()->SetRangeUser(0,1000);
    hWishboneEProj[true]->Draw();
    hWishboneEProj[false]->Draw("Same");
    drawHLine(0., defaultCanvas, 1);
    printCanvas("WishboneEnergy");
    
    hWishboneTProj->SetMinimum(0);
    hWishboneTProj->SetMaximum(20);
    hWishboneTProj->Draw("E0");
    drawVLine(T_p_lo/1000., defaultCanvas, 2);
    drawVLine(T_p_hi/1000., defaultCanvas, 2);
    drawVLine(T_p_min/1000., defaultCanvas, 4);
    drawVLine(T_p_max/1000., defaultCanvas, 4);
    printCanvas("WishboneTime");
    
    hWishboneFiducialTProj->Draw();
    printCanvas("WishboneFiducialTime");
    
    if(isCombined) {
        GausWishboneFit GWF("WishboneFit", myA);
        GWF.setWishbone(hWishboneBGSub);
        GWF.fitModel();
        GWF.extractAsymmetry();
        GWF.calcGapFill();
        
        WishboneMidcentroid WM("WishboneMidcentroid", myA);
        WM.setWishbone(hWishboneBGSub);
        WM.extractAsymmetry();
        GWF.compareAsym(WM);
    }
}
