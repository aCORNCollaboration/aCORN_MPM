/// \file PulserPlugin.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2016

#include "PulserPlugin.hh"
#include "GraphUtils.hh"
#include "GraphicsUtils.hh"

PulserPlugin::PulserPlugin(RunAccumulator* RA, const string& nm, const string& inflname):
RunAccumulatorPlugin(RA, nm, inflname),
hPulserSignal(this) {
    
    setAnalysisCuts();
    ignoreMissingHistos = true;
    
    hPulserTiming = registerSavedHist("hPulserTiming", "proton pulser timing correlation", 400, -10, 10);
    hPulserTiming->GetXaxis()->SetTitle("(time to preceding pulses) - 1s [ms]");
    
    hPulserTimingPrecis = registerSavedHist("hPulserTimingPrecis", "proton pulser timing correlation", 200, -100, 100);
    hPulserTimingPrecis->GetXaxis()->SetTitle("(time to preceding pulses) - t_{0} [#mus]");
    
    TH1F hPulserSignalTemplate("hPulserSignal", "proton pulser signals", 200, 0, 20);
    hPulserSignalTemplate.GetXaxis()->SetTitle("proton detector ADC channels (#times 10^{3})");
    initRegions(hPulserSignal);
    hPulserSignal.setTemplate(hPulserSignalTemplate);
    
    TH1* hPTimingGapTemplate = logHist("hPTimingGapTemplate", "proton events timing gap", 112, 1e-6, 100);
    hPTimingGap = registerSavedHist("hPTimingGap", *hPTimingGapTemplate);
    hPTimingGap->GetXaxis()->SetTitle("gap between proton events [s]");
    delete hPTimingGapTemplate;
    
    vGapTime = registerNamedVector("vGapTime", 1);
}

void PulserPlugin::setAnalysisCuts() {
    DataMode dm = myA->dataMode;
    if(dm == NGC) {
        E_puls_lo = 3000;
        E_puls_hi = 6500;
        T_pulser = 0.99580e9;
    } else if(dm == NG6) {
        E_puls_lo = 8000;
        E_puls_hi = 12000;
        T_pulser = 0.99585e9;
    }
}
   
void PulserPlugin::initRegions(FGBGRegionsHist& h) {
    h.addRegion(-1e5, +1e5, true);
    h.addRegion(-2e6, -1e6, false);
}

void PulserPlugin::fillCoreHists(BaseDataScanner& PDS, double weight) {

    // consider only first record for each proton
    if(PDS.E_p <= 0) return;
    
    // in pulser energy range?
    bool isPulserEnergy = E_puls_lo < PDS.E_p && PDS.E_p < E_puls_hi;
    // truncate history of previous pulses to up to 1.011s ago (slightly more than "coarse" search window)
    while(prevPulses.size() && (prevPulses.front() > PDS.T_p || prevPulses.front() < PDS.T_p - 1.011e9)) prevPulses.pop_front();
        
    // timing coincidences with previous pulses
    for(auto pp: prevPulses) {
        double dt = PDS.T_p-pp-1e9;
        if(dt < -1e7) break; // stop evaluating outside coarse search range
        double dtp = PDS.T_p-pp-T_pulser;
        hPulserSignal.fill(dtp, PDS.E_p_0/1000., weight);
        if(isPulserEnergy) {
            hPulserTiming->Fill(dt*1e-6, weight);
            hPulserTimingPrecis->Fill(dtp*1e-3,weight);
        }
    }
        
    // record candidate pulse to previous pulses list
    if(isPulserEnergy) prevPulses.push_back(PDS.T_p);
    
    // proton timing gap
    double dgap = (PDS.T_p - T_prev_p)*1e-9;
    if(dgap > 0) hPTimingGap->Fill(dgap, weight);
    if(dgap > tau_long) (*vGapTime)[0] += dgap;
    T_prev_p = PDS.T_p;
}

void PulserPlugin::calculateResults() {
    hPulserSignal.makeRates(1);
    isCalculated = true;
}

double log_binsize(const TH1& h) {
    const TAxis* Ax = h.GetXaxis();
    int nb = Ax->GetNbins();
    return log(Ax->GetBinUpEdge(nb)/Ax->GetBinLowEdge(1))/nb;
}

void PulserPlugin::makeAnaResults() {
    AnaResult baseResult = myA->makeBaseResult();

    // coarse timing
    double pt0 = hPulserTiming->GetBinCenter(hPulserTiming->GetMaximumBin()); // [ms]
    double offset_from_expected = pt0 + 1e3 - T_pulser*1e-6; // [ms]
    printf("Coarse timing offset by %g ms from expected.\n", offset_from_expected);
    
    // fine timing
    if(fabs(offset_from_expected) < 0.1) {
        TF1* gbgfit = new TF1("gbgfit", "[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2])) + [3]", -100, 100);
        addDeletable(gbgfit);
        gbgfit->SetParameter(0, hPulserTimingPrecis->GetMaximum());
        gbgfit->SetParameter(1, hPulserTimingPrecis->GetBinCenter(hPulserTimingPrecis->GetMaximumBin()));
        gbgfit->SetParameter(2, 20.);
        hPulserTimingPrecis->Fit(gbgfit,"R");
        
        baseResult.value = gbgfit->GetParameter(2)/1000.;
        baseResult.err = gbgfit->GetParError(2)/1000.;
        AcornDB::ADB().uploadAnaResult("pulser_jitter", "pulser timing jitter [ms]", baseResult);
        
        baseResult.value = T_pulser*1e-9 + gbgfit->GetParameter(1)*1e-6;
        baseResult.err = gbgfit->GetParError(1)*1e-6;
    } else {
        baseResult.value = pt0*1e-3 + 1.0;
        baseResult.err = 0;
    }
    AcornDB::ADB().uploadAnaResult("pulser_timing", "pulser timing correlation peak [s]", baseResult);
    
    TF1* pulpkft = new TF1("pulpkft", "gaus", E_puls_lo/1000, E_puls_hi/1000);
    addDeletable(pulpkft);
    hPulserSignal.hRates[true]->Fit(pulpkft,"R");
    baseResult.value = pulpkft->GetParameter(1)*1000.;
    baseResult.err = pulpkft->GetParError(1)*1000.;
    AcornDB::ADB().uploadAnaResult("pulser_center", "pulser signal center [channels]", baseResult);
    baseResult.value = pulpkft->GetParameter(2)*1000.;
    baseResult.err = pulpkft->GetParError(2)*1000.;
    AcornDB::ADB().uploadAnaResult("pulser_width", "pulser signal width [channels]", baseResult);
    
    baseResult.value = integralAndErrorInterp(hPulserSignal.hRates[true], E_puls_lo/1000, E_puls_hi/1000, baseResult.err, true);
    AcornDB::ADB().uploadAnaResult("pulser_rate", "pulser paired-event rate [Hz]", baseResult);
    
    
    // calculation of efficiency for catching random proton events
    
    TF1* pdtfit = new TF1("pdtfit", "[0]*x*exp(-x/[1])", 5e-4, 0.05);
    addDeletable(pdtfit);
    double tau0 = hPTimingGap->GetBinCenter(hPTimingGap->GetMaximumBin());
    if(tau0 < 1e-3 || tau0 > 1e-1) tau0 = 0.01;
    pdtfit->SetParameter(1,tau0);
    pdtfit->SetParameter(0,1);
    double h =  hPTimingGap->GetMaximum()/pdtfit->Eval(tau0);
    pdtfit->SetParameter(0, h);
    printf("Initial guess: tau0 = %g, h = %g\n", tau0, h);
    hPTimingGap->Fit(pdtfit, "R");
    
    h = pdtfit->GetParameter(0);
    
    baseResult.value = tau0 = pdtfit->GetParameter(1);
    baseResult.err = pdtfit->GetParError(1);
    AcornDB::ADB().uploadAnaResult("proton_random_dt", "uncorrelated proton delta time [s]", baseResult);
    
    double fitcounts = h * tau0 / log_binsize(*hPTimingGap);
    double expcounts = myA->runTimes.total()/tau0;
    baseResult.value = fitcounts/expcounts;
    baseResult.err = 0;
    AcornDB::ADB().uploadAnaResult("proton_random_effic", "uncorrelated proton DAQ efficiency", baseResult);
    
    // dead fraction from summing proton gaps >= 50ms, less distribution-expected value
    baseResult.value = ((*vGapTime)[0] - fitcounts*(tau_long+tau0)*exp(-tau_long/tau0)) / myA->runTimes.total();
    baseResult.err = 0;
    AcornDB::ADB().uploadAnaResult("proton_gap_loss", "DAQ dead fraction from proton gaps", baseResult);
    // number of gap events, less expected
    baseResult.value = hPTimingGap->Integral(hPTimingGap->FindBin(tau_long), hPTimingGap->GetNbinsX()+1) - fitcounts*exp(-tau_long/tau0);
    baseResult.err = sqrt(baseResult.value);
    AcornDB::ADB().uploadAnaResult("n_gaps", "number of gaps in proton data", baseResult);
}

void PulserPlugin::makePlots() {
    assert(isCalculated);
    defaultCanvas->SetRightMargin(0.04);
    defaultCanvas->SetLogy(true);

    hPulserSignal.hRates[false]->SetMinimum(1e-3);
    hPulserSignal.hRates[false]->SetMaximum(2);
    hPulserSignal.hRates[false]->Draw();
    hPulserSignal.hRates[true]->Draw("Same");
    addDeletable(drawVLine(E_puls_lo/1000., defaultCanvas, 1));
    addDeletable(drawVLine(E_puls_hi/1000., defaultCanvas, 1));
    printCanvas("PulserSignal");
  
    hPulserTiming->Draw();
    addDeletable(drawVLine((T_pulser-1e9)*1e-6, defaultCanvas, 2));
    printCanvas("PulserTiming");
    hPulserTimingPrecis->Draw();
    printCanvas("PulserTimingPrecision");
    
    defaultCanvas->SetLogx(true);
    hPTimingGap->GetXaxis()->SetRangeUser(1e-6,1);
    hPTimingGap->Draw();
    printCanvas("ProtonTimingGap");
    
    defaultCanvas->SetLogx(false);
    defaultCanvas->SetLogy(false);
}
