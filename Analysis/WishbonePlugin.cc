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
#include "IterRangeFitter.hh"

#include <TStyle.h>
#include <TColor.h>
#include <cassert>

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
    
    setAnalysisCuts();
    
    wbTimingRegions.newRegion();
    wbTimingRegions.appendRange(T_p_min, T_p_lo);
    wbTimingRegions.appendRange(T_p_hi, T_p_max);
    wbTimingRegions.newRegion();
    wbTimingRegions.appendRange(T_p_lo, T_p_hi);
    wbTimingRegions.scaleNormsTo(1);
    
    TH1F hProtonTemplate("hProtonSignal", "Proton detector signal", 200,0,20);
    hProtonTemplate.GetXaxis()->SetTitle("proton detector ADC channels (#times 10^{3})");
    hProtonSignal.setTemplate(hProtonTemplate, wbTimingRegions);
  
    TH1F hNVetoTemplate("hNVeto","Veto panels count",9,-0.5,8.5);
    hNVetoTemplate.GetXaxis()->SetTitle("Number of veto panels");
    hNVeto.setTemplate(hNVetoTemplate, wbTimingRegions);
    
    TH1F hVetoSumTemplate("hVetoSum","Veto PMTs signal sum",100,0,20);
    hVetoSumTemplate.GetXaxis()->SetTitle("raw signal sum (#times 10^{3})");
    hVetoSum.setTemplate(hVetoSumTemplate, wbTimingRegions);
    
    unsigned int nEnBins = 400;
    double E0 = 0;
    double E1 = 1600;
    double timeBin = 0.02; // time bin width, microseconds
    
    // note 5ns shift in time axis bins, to avoid bin boundaries at 10ns+/-error unit boundaries.
    TH2F hTemplate("hTemplate","aCORN Wishbone", nEnBins, E0, E1, 10/timeBin + 1, - 0.005, 10. + timeBin - 0.005);
    hTemplate.GetXaxis()->SetTitle("Electron energy [keV]");
    hTemplate.GetYaxis()->SetTitle("Proton TOF [#mus]");
    hWishbone = (TH2F*)registerSavedHist("hWishbone",hTemplate);
    //if(myA->dataMode == NG6) hWishbone->SetTitle("aCORN NG-6 Wishbone");
    if(myA->dataMode == NGC) hWishbone->SetTitle("aCORN NG-C Wishbone");
    hWishbone->GetYaxis()->SetRange();
    
    TH2F hNETemplate("hNE","PMT trigger counts", 100, 0, 1000, 20, -0.5, 19.5);
    hNETemplate.GetYaxis()->SetTitle("Number of PMTs triggered");
    hNETemplate.GetXaxis()->SetTitle("Electron energy [keV]");
    hNE.setTemplate(hNETemplate, wbTimingRegions);
    
    TH2F hPSTemplate("hChanSpec","PMT Spectra", 200, 0, (1<<15)/1000., NCH_MAX, -0.5, NCH_MAX-0.5);
    hPSTemplate.GetXaxis()->SetTitle("raw signal (#times 10^{3})");
    hPSTemplate.GetYaxis()->SetTitle("detector channel number");
    hChanSpec.setTemplate(hPSTemplate, wbTimingRegions);
    
    TH2F hMultTemplate("hModuleMult","Module Multiplicity", CHAN_PER_MOD+1, -0.5, CHAN_PER_MOD+0.5, CHAN_PER_MOD+1, -0.5, CHAN_PER_MOD+0.5);
    hMultTemplate.GetXaxis()->SetTitle("Module 1 Triggers");
    hMultTemplate.GetYaxis()->SetTitle("Module 2 Triggers");
    hModuleMult.setTemplate(hMultTemplate, wbTimingRegions);
    
    TH2F hPosTemplate("hPos","Beta hit positions", 100, -3, 3, 100, -3, 3);
    hPosTemplate.GetXaxis()->SetTitle("x [PMT spacings]");
    hPosTemplate.GetYaxis()->SetTitle("y [PMT spacings]");
    hPos.setTemplate(hPosTemplate, wbTimingRegions);
    
    TH2F hPosSigmaTemplate("hPosSigma","Beta hit position spread", 100, 0, 2, 100, 0, 2);
    hPosSigmaTemplate.GetXaxis()->SetTitle("#sigma_{x} [PMT spacings]");
    hPosSigmaTemplate.GetYaxis()->SetTitle("#sigma_{y} [PMT spacings]");
    hPosSigma.setTemplate(hPosSigmaTemplate, wbTimingRegions);
    
    
    TH2F hEnergyRadiusTemplate("hEnergyRadius","Beta hit positions", 100, 0, 4, 100, 0, 800);
    hEnergyRadiusTemplate.GetXaxis()->SetTitle("radius^{2} [(PMT spacings)^{2}]");
    hEnergyRadiusTemplate.GetYaxis()->SetTitle("Energy [keV]");
    hEnergyRadius.setTemplate(hEnergyRadiusTemplate, wbTimingRegions);
    
    ignoreMissingHistos = true;
    TRateCategories hRateHistoryTemplate("hRateHistory","Event rate",1./6.);
    hRateHistory = (TRateCategories*)registerObject("hRateHistory", hRateHistoryTemplate);
}

void WishbonePlugin::setAnalysisCuts() {
    DataMode dm = myA->dataMode;
    if(dm == NGC) {
        E_p_lo = 550;
        E_p_hi = 1500;
        T_p_min = 750;
        T_p_lo = 3000;
        T_p_hi = 4500;
        T_p_max = 9500;
    } else if(dm == NG6) {
        E_p_lo = 650;
        E_p_hi = 2400;
        T_p_min = 750;
        T_p_lo = 2750;
        T_p_hi = 4500;
        T_p_max = 9500;
    }
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
    Double_t t_proton_h = (PDS.T_0+PDS.T_p)*1e-9/3600.;
    hRateHistory->AddPoint(RATE_PROTON, t_proton_h, weight);
    
    // beta decay protons
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
        
        if(200 < PDS.E_recon && PDS.E_recon < 600) {
            for(int fg = 0; fg < 2; fg ++) {
                for(auto cut: wbTimingRegions.regions[fg]) { 
                    if(cut.first < PDS.T_e2p && PDS.T_e2p <= cut.second) {
                        hRateHistory->AddPoint(fg? RATE_WB_FG : RATE_WB_BG, t_proton_h, weight);
                        if(PDS.nV) hRateHistory->AddPoint(fg? RATE_VETO_FG : RATE_VETO_BG, t_proton_h, weight);
                    }
                }
            }
        }
    }
}

void WishbonePlugin::calculateResults() {
    hWishboneEProj[false] = hWishboneEProj[true] = NULL;
    TH2Slicer wbs(hWishbone);
    double tfg = wbs.projSlice(T_p_lo/1000., T_p_hi/1000., hWishboneEProj[true]);
    double tbg = wbs.projSlice(T_p_min/1000., T_p_lo/1000., hWishboneEProj[false]);
    tbg += wbs.projSlice(T_p_hi/1000., T_p_max/1000., hWishboneEProj[false]);
    
    hWishboneBGSub = wbs.subtractProfile(hWishboneEProj[false],1./tbg);
    
    double s0 = 1000./myA->runTimes->GetTotal()/hWishboneEProj[true]->GetBinWidth(1);
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
    normalize_to_bin_width(hWishboneTProj, 1./myA->runTimes->GetTotal());
    hWishboneTProj->GetYaxis()->SetTitle("event rate [Hz/#mus]");
    
    hWishboneBGSub->Scale(s0/hWishboneTProj->GetBinWidth(1));
    hWishboneBGSub->GetZaxis()->SetTitle("rate [Hz/MeV/#mus]");
    
    int eb0 = hWishboneBGSub->GetXaxis()->FindBin(100);
    int eb1 = hWishboneBGSub->GetXaxis()->FindBin(300);
    hWishboneFiducialTProj = hWishboneBGSub->ProjectionY("_tproj", eb0, eb1, "e");
    hWishboneFiducialTProj->Scale(hWishboneBGSub->GetXaxis()->GetBinWidth(1)/1000.);
    hWishboneFiducialTProj->GetYaxis()->SetTitle("event rate [Hz/#mus]");

    hWBRate = (TH2F*)hWishbone->Clone("hWBRate");
    normalize_to_bin_area(hWBRate, 1000./myA->runTimes->GetTotal());
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
    
    hRateHistory->SummarizeWindow();
    
    isCalculated = true;
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
    
    // Energy projection endpoint fit and rescaled energy range rate integrals
    IterRangeErfc EF(650,75);
    EF.nsigmalo = 3;
    EF.doFit(hWishboneEProj[true]);
    baseResult.value = EF.myF->GetParameter(1);
    baseResult.err = EF.myF->GetParError(1);
    AcornDB::ADB().uploadAnaResult("spectrum_edge_center", "fit to center of spectrum upper edge", baseResult);
    baseResult.value = EF.myF->GetParameter(2);
    baseResult.err = EF.myF->GetParError(2);
    AcornDB::ADB().uploadAnaResult("spectrum_edge_width", "fit to width of spectrum upper edge", baseResult);
    Escale = 650./EF.myF->GetParameter(1); // energy re-scaling factor
    double rE0 = 100/Escale;
    double rE1 = 900/Escale;
    baseResult.value = integralAndErrorInterp(hWishboneEProj[true], rE0, rE1, baseResult.err, true) / 1000.;
    baseResult.err /= 1000;
    AcornDB::ADB().uploadAnaResult("wb_erange_rate", "wishbone rate in ~100--900keV [Hz]", baseResult);
    baseResult.value = integralAndErrorInterp(hWishboneEProj[false], rE0, rE1, baseResult.err, true) / 1000.;
    baseResult.err /= 1000;
    AcornDB::ADB().uploadAnaResult("bg_erange_rate", "background rate in ~100--900keV [Hz]", baseResult);
    
    // raw proton signal peak
    TF1* fPPkFit = new TF1("fPPkFit", "[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2])) + (1+[3]*x)/([6] + [4]*x + [5]*x*x+[7]*x*x*x)", E_p_lo/1000,  0.5+E_p_hi/1000);
    double c0 = 0.959; //0.5*(E_p_lo + E_p_hi)/1000.;
    TH1* hPPk = hProtonSignal.hRates[false];
    double ph = hPPk->GetBinContent(hPPk->FindBin(c0));
    fPPkFit->SetParameter(0,ph);
    fPPkFit->FixParameter(1,c0);
    fPPkFit->FixParameter(2,0.123);
    for(int i=3; i<=6; i++) fPPkFit->SetParameter(i,0);
    fPPkFit->SetParameter(6,3/ph);

    hPPk->Fit(fPPkFit,"R");
    ph = fPPkFit->GetParameter(0);
    double ps = fPPkFit->GetParameter(2);
    printf("Proton peak rate: %g Hz\n", ph*ps*sqrt(2*M_PI));
    
    ManualWishboneSeparator MWS("WBFit", myA);
    MWS.setWishbone(hWishboneBGSub);
    MWS.extractAsymmetry();
    baseResult.value = MWS.fitAsym(50,400,baseResult.err);
    AcornDB::ADB().uploadAnaResult("obs_asym", "observed uninterpreted asymmetry", baseResult);
}

void WishbonePlugin::makePlots() {
    assert(isCalculated);
    defaultCanvas.SetRightMargin(0.20);

    bool isCombined = myA->runTimes->Size() > 1000;
    if(true || !isCombined) {
        hWishboneBGSub->Rebin2D(2,2);
        hWishboneBGSub->Scale(1./4.);
    }
    hWishboneBGSub->SetMinimum(-10);
    hWishboneBGSub->SetMaximum(10);
    hWishboneBGSub->GetXaxis()->SetRangeUser(0,1000);
    hWishboneBGSub->GetYaxis()->SetRangeUser(2,5);
    makeRBpalette();
    hWishboneBGSub->Draw("Col Z");
    addDeletable(drawHLine(T_p_lo/1000., &defaultCanvas, 2));
    addDeletable(drawHLine(T_p_hi/1000., &defaultCanvas, 2));
    addDeletable(drawHLine(T_p_min/1000., &defaultCanvas, 4));
    addDeletable(drawHLine(T_p_max/1000., &defaultCanvas, 4));
    printCanvas("Wishbone");
    
    defaultCanvas.SetRightMargin(0.15);
    
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
    
    defaultCanvas.SetLogz(true);
    
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
    
    defaultCanvas.SetLogz(false);
    
    defaultCanvas.SetRightMargin(0.04);
    
    defaultCanvas.SetLogy(true);

    hNVeto.hRates[true]->SetMinimum(1e-7);
    hNVeto.hRates[true]->SetMaximum(10);
    hNVeto.hRates[true]->Draw();
    hNVeto.hRates[false]->Draw("Same");
    printCanvas("NVeto");
    
    hProtonSignal.hRates[false]->SetMinimum(isCombined?1e-4:1e-3);
    hProtonSignal.hRates[false]->SetMaximum(50);
    hProtonSignal.hRates[false]->Draw();
    hProtonSignal.hRates[true]->Draw("Same");
    addDeletable(drawVLine(E_p_lo/1000., &defaultCanvas, 2));
    addDeletable(drawVLine(E_p_hi/1000., &defaultCanvas, 2));
    printCanvas("ProtonSignal");
      
    defaultCanvas.SetLogy(false);
    
    hVetoSum.hRates[false]->SetMinimum(-20);
    hVetoSum.hRates[false]->SetMaximum(200);
    hVetoSum.hRates[false]->Draw();
    addDeletable(drawHLine(0, &defaultCanvas, 1));
    hVetoSum.hRates[true]->Draw("Same");
    printCanvas("VetoSum");
    
    hWishboneEProj[true]->SetMinimum(myA->dataMode == NG6? -0.2 : -1);
    hWishboneEProj[true]->SetMaximum(myA->dataMode == NG6? 3 : 15);
    hWishboneEProj[true]->GetXaxis()->SetRangeUser(0,1000);
    hWishboneEProj[true]->Draw();
    hWishboneEProj[false]->Draw("Same");
    addDeletable(drawVLine(100/Escale, &defaultCanvas, 1, 2));
    addDeletable(drawVLine(900/Escale, &defaultCanvas, 1, 2));
    addDeletable(drawHLine(0., &defaultCanvas, 1));
    printCanvas("WishboneEnergy");
    
    hWishboneTProj->SetMinimum(0);
    hWishboneTProj->SetMaximum(20);
    hWishboneTProj->Draw("E0");
    addDeletable(drawVLine(T_p_lo/1000., &defaultCanvas, 2));
    addDeletable(drawVLine(T_p_hi/1000., &defaultCanvas, 2));
    addDeletable(drawVLine(T_p_min/1000., &defaultCanvas, 4));
    addDeletable(drawVLine(T_p_max/1000., &defaultCanvas, 4));
    printCanvas("WishboneTime");
   
    // scaled background and background-subtracted rate history plots
    assert(hRateHistory);
    vector< pair<Int_t, Double_t> > rs;
    rs.push_back(pair<Int_t, Double_t>(RATE_WB_BG, 1./3600/wbTimingRegions.norms[0]));
    TGraph* g0 = hRateHistory->MakeGraph(rs);
    rs[0].first = RATE_VETO_BG;
    TGraph* g2 = hRateHistory->MakeGraph(rs);
    
    rs[0].second *= -1.0;
    rs.push_back(pair<Int_t, Double_t>(RATE_VETO_FG, 1./3600));
    TGraph* g3 = hRateHistory->MakeGraph(rs);
    rs[0].first = RATE_WB_BG;
    rs[1].first = RATE_WB_FG;
    TGraph* g1 = hRateHistory->MakeGraph(rs);
    
    g1->SetLineColor(2);
    g0->SetLineColor(4);
    
    g0->SetMinimum(-0.1);
    g0->SetTitle("Event rate");
    g0->Draw("ALZ");
    g0->GetXaxis()->SetTitle("run time [h]");
    g0->GetYaxis()->SetTitle("rate [Hz]");
    if(g0->GetN()) {
        double x0 = g0->GetX()[0];
        g0->GetXaxis()->SetRangeUser(x0 > 1? x0 : 0, g0->GetX()[g0->GetN()-1]);
    }
    g1->Draw("LZ");
    
    g2->SetLineColor(7);
    g2->Draw("LZ");
    g3->SetLineColor(6);
    g3->Draw("LZ");
    
    printCanvas("RateHistory");
    delete g0;
    delete g1;
    delete g2;
    delete g3;
    
    /////////////////////////////////////////////////////////////////
    
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
