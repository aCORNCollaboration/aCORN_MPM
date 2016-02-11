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

    if(PDS.E_p_0 <= 0) return;
    
    // proton detector pulser signal
    // note: using E_p instead of E_p_0 to only get first occurrence of each proton
    if(PDS.E_p > 0) {
        bool isPulserEnergy = E_puls_lo < PDS.E_p && PDS.E_p < E_puls_hi; // in pulser energy range?
        if(isPulserEnergy) {
            prevPulses.push_back(PDS.T_p);
            while(prevPulses.front() > PDS.T_p || prevPulses.front() < PDS.T_p - 1.2e9) prevPulses.pop_front();
        }
        for(auto pp: prevPulses) {
            double dt = PDS.T_p-pp-1e9;
            double dtp = PDS.T_p-pp-T_pulser;
            if(isPulserEnergy) {
                hPulserTiming->Fill(dt*1e-6, weight);
                hPulserTimingPrecis->Fill(dtp*1e-3,weight);
            }
            if(dt < -1e7) break;
            hPulserSignal.fill(dtp, PDS.E_p_0/1000., weight);
        }
    }
}

void PulserPlugin::calculateResults() {
    hPulserSignal.makeRates(1);
    isCalculated = true;
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
        gbgfit->SetParameter(0, hPulserTimingPrecis->GetMaximum());
        gbgfit->SetParameter(1, hPulserTimingPrecis->GetBinCenter(hPulserTimingPrecis->GetMaximumBin()));
        gbgfit->SetParameter(2, 20.);
        hPulserTimingPrecis->Fit(gbgfit,"R+");
        baseResult.value = T_pulser*1e-9 + gbgfit->GetParameter(1)*1e-6;
        baseResult.err = gbgfit->GetParError(1)*1e-6;
    } else {
        baseResult.value = pt0*1e-3 + 1.0;
        baseResult.err = 0;
    }
    AcornDB::ADB().uploadAnaResult("pulser_timing", "pulser timing correlation peak [s]", baseResult);
    
    baseResult.value = integralAndErrorInterp(hPulserSignal.hRates[true], E_puls_lo/1000, E_puls_hi/1000, baseResult.err, true);
    AcornDB::ADB().uploadAnaResult("pulser_rate", "pulser paired-event rate [Hz]", baseResult);
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
    
    defaultCanvas->SetLogy(false);
}
