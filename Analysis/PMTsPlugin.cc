/// \file PMTsPlugin.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "PMTsPlugin.hh"
#include <TStyle.h>

PMTsPlugin::PMTsPlugin(RunAccumulator* RA, OutputManager* pnt, const string& nm, const string& inflname):
RunAccumulatorPlugin(RA, pnt, nm, inflname) {
    unsigned int nEnBins = 400;
    double E0 = 0;
    double E1 = 3000;
    
    
    hEnergy = registerSavedHist("hEnergy","Electron energy spectrum", nEnBins, E0, E1);
    hEnergy->GetXaxis()->SetTitle("Electron energy [keV]");
    
    TH2F hNETemplate("hNEtemplate","PMT trigger counts", 100, 0, 3000, 20, -0.5, 19.5);
    hNETemplate.GetYaxis()->SetTitle("Number of PMTs triggered");
    hNETemplate.GetXaxis()->SetTitle("Electron energy [keV]");
    hNE = (TH2*)registerSavedHist("hNE",hNETemplate);
    
    hEETime =  registerSavedHist("hEETime","Electron event timing difference", 999, 0, 10);
    hEETime->GetXaxis()->SetTitle("time difference [#mus]");
    
    TH2F hPSTemplate("hChanSpecT","PMT Spectra", 200, 0, 20., NCH_MAX, -0.5, NCH_MAX-0.5);
    hPSTemplate.GetXaxis()->SetTitle("raw signal (#times 10^{3})");
    hPSTemplate.GetYaxis()->SetTitle("detector channel number");
    hPSTemplate.GetYaxis()->SetTitleOffset(1.4);
    hChanSpec = (TH2*)registerSavedHist("hChanSpec",hPSTemplate);
}


void PMTsPlugin::fillCoreHists(BaseDataScanner& PDS, double weight) {
    double vetosum = 0;
    for(unsigned int i=0; i<NCH_MAX; i++) {
        if(PDS.E_PMT[i]) {
            hChanSpec->Fill(PDS.E_PMT[i]/1000., i, weight);
            if(i<N_V_PMT) vetosum += PDS.E_PMT[i];
        }
    }
    
    //if(PDS.E_p == 0) hEETime->Fill((PDS.T_p - prev_e_time)/1000);
    //prev_e_time = PDS.T_p;
    
    if(PDS.E_p <= 0) hEETime->Fill( (PDS.T_e2p - prev_e_time)/1000 );
    prev_e_time = PDS.T_e2p;
    
    //if(vetosum) hVetoSum.fill(PDS.T_e2p, vetosum/1000., weight);
    
    hEnergy->Fill(PDS.E_recon, weight);
    hNE->Fill(PDS.E_recon, PDS.nE, weight);
}

void PMTsPlugin::makePlots() {
    gStyle->SetPalette(1);
    
    hEETime->Draw("HIST");
    gPad->SetLogy(true);
    myA->printCanvas("hEETime");
    
    TH2* hNERate = (TH2*)myA->hToRate(hNE,0);
    hNERate->SetMinimum(0.01);
    hNERate->Draw("Col Z");
    gPad->SetLogz(true);
    myA->printCanvas("hNE");
    
    TH1* hEnergyR = myA->hToRate(hEnergy,1);
    hEnergyR->GetYaxis()->SetTitle("event rate [Hz/keV]");
    hEnergyR->SetMaximum(2*hEnergyR->GetBinContent(hEnergyR->FindBin(200)));
    hEnergyR->Draw();
    myA->printCanvas("hEnergy");
}
