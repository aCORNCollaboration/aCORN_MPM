#include "SourceCalPlugin.hh"
#include "AcornDB.hh"
#include <TMath.h>
#include <TF1.h>
#include <TLatex.h>

SourceCalPlugin::SourceCalPlugin(RunAccumulator* RA): AnalyzerPlugin(RA,"SourceCal") {    
    hEnergy = registerHist("hEnergy", "Calibrated energy", 200, 0, 1500);
    hEnergy->GetXaxis()->SetTitle("Energy [keV]");
    
    hEnergyRecal = registerHist("hEnergyRecal", "Re-calibrated energy", 200, 0, 1500);
    hEnergyRecal->GetXaxis()->SetTitle("Energy [keV]");
    
    for(unsigned int i=0; i<N_E_PMT; i++) {
    //for(unsigned int i=9; i<11; i++) {
        hPMTSig[i] = registerHist("hPMTSig_"+itos(i), "PMT " + itos(i+1) + " signal", 200, 0, 15);
        hPMTSig[i]->GetXaxis()->SetTitle("ADC channels [#times 10^{3}]");
    }
}

void SourceCalPlugin::fillCoreHists(BaseDataScanner& PDS, double weight) {
    if(PDS.Pos.r2() > 0.5*0.5) return;
    hEnergy->Fill(PDS.E_recon, weight);
    hEnergyRecal->Fill(PDS.getCal()->calEnergy(PDS.E_PMT + N_V_PMT),weight);
    if(850 < PDS.E_recon && PDS.E_recon < 1200) {
        //for(unsigned int i=9; i<11; i++)
        for(unsigned int i=0; i<N_E_PMT; i++)
            hPMTSig[i]->Fill(PDS.E_PMT[i+N_V_PMT]/1000., weight);
    }
}

Double_t poissonf(Double_t* x, Double_t* par) {
    return par[0]*TMath::Poisson(x[0]*par[1], par[2]);
}

void fit_pks(TH1* hPk, double& p1, double& p2, bool show_fit = true) {
    TF1 fGaus("fGaus","gaus",0,12);
    TF1 fPois("pois", &poissonf, 0, 12, 3);
    fPois.SetLineColor(1);
        
    double c0 = hPk->GetBinCenter(hPk->GetMaximumBin());
    fGaus.SetRange(c0-0.3*c0,c0+0.3*c0);
    hPk->Fit(&fGaus,"QRN");
    
    c0 = fGaus.GetParameter(1);
    double w = fGaus.GetParameter(2);
    double N = c0*c0/(w*w);
        
    std::cout << "Estimating N = " << N << " from width w = " << w << " out of c = " << c0 << "\n";
    fPois.SetParameter(0, 1.0);
    fPois.SetParameter(1, N/c0);
    fPois.SetParameter(2, N);
    std::cout << "Initial height = " << fPois.Eval(c0) << "\n";
    fPois.SetParameter(0,fGaus.GetParameter(0)/fPois.Eval(c0));
    fPois.SetRange(c0-1.5*w,c0+1.5*w);
        
    hPk->Fit(&fPois,show_fit?"R":"QRN");
    
    p1 = fPois.GetParameter(1);
    p2 = fPois.GetParameter(2);
}

void SourceCalPlugin::bgSubtrPlots(SourceCalPlugin& bg) {
    TH1* hEnergyRate = myA->hToRate(hEnergy,1);
    hEnergyRate->Scale(1000);
    hEnergyRate->GetYaxis()->SetTitle("event rate [mHz/keV]");
    hEnergyRate->GetYaxis()->SetTitleOffset(1.4);
    hEnergyRate->SetLineColor(2);
    TH1* hEnergyRateBG = bg.myA->hToRate(bg.hEnergy,1);
    hEnergyRateBG->Scale(1000);
    hEnergyRateBG->SetLineColor(4);
    
    myA->defaultCanvas->SetLogy(true);
    hEnergyRate->Draw();
    hEnergyRateBG->Draw("Same");
    myA->printCanvas("Energy_BG");
    
    myA->defaultCanvas->SetLogy(false);
    hEnergyRate->Add(hEnergyRateBG,-1.0);
    hEnergyRate->SetMaximum(80);
    
    TF1 fPois("pois", &poissonf, 900, 1150, 3);
    fPois.SetParameter(0,60);
    fPois.SetParameter(0,1000/257.);
    fPois.SetParameter(2,257.);
    fPois.SetLineColor(1);
    hEnergyRate->Fit(&fPois,"R");
    printf("Original energy: %g PE/MeV\n",fPois.GetParameter(2)/0.9948);
    
    hEnergyRate->Draw();
    hEnergyRateBG->Draw("Same");
    myA->printCanvas("Energy");
    
    //------------------------
    
    TH1* hEnergyRecalRate = myA->hToRate(hEnergyRecal,1);
    hEnergyRecalRate->Scale(1000);
    hEnergyRecalRate->GetYaxis()->SetTitle("event rate [mHz/keV]");
    hEnergyRecalRate->GetYaxis()->SetTitleOffset(1.4);
    hEnergyRecalRate->SetLineColor(2);
    TH1* hEnergyRecalRateBG = bg.myA->hToRate(bg.hEnergyRecal,1);
    hEnergyRecalRateBG->Scale(1000);
    hEnergyRecalRateBG->SetLineColor(4);
    
    hEnergyRecalRate->Add(hEnergyRecalRateBG,-1.0);
    hEnergyRecalRate->SetMaximum(80);
    TF1 fPoisR("pois", &poissonf, 1000, 1200, 3);
    fPoisR.SetParameter(0,60);
    fPoisR.SetParameter(0,1000/257.);
    fPoisR.SetParameter(2,257.);
    fPoisR.SetLineColor(1);
    hEnergyRecalRate->Fit(&fPoisR,"R");
    printf("Recalibrated energy: %g PE/MeV\n",fPoisR.GetParameter(2)/0.9948);
    
    hEnergyRecalRate->Draw();
    hEnergyRecalRateBG->Draw("Same");
    myA->printCanvas("EnergyRecal");
    
    //---------------------
    
    double pmt_res[N_E_PMT];
    double pmt_n0[N_E_PMT];
    //for(unsigned int i=9; i<11; i++) {
    for(unsigned int i=0; i<N_E_PMT; i++) {
        TH1* hPMTRate = myA->hToRate(hPMTSig[i],1);
        hPMTRate->GetYaxis()->SetTitle("event rate [mHz/channel]");
        hPMTRate->GetYaxis()->SetTitleOffset(1.4);
        hPMTRate->SetLineColor(2);
        TH1* hPMTRateBG = bg.myA->hToRate(bg.hPMTSig[i],1);
        hPMTRateBG->SetLineColor(4);
        
        hPMTRate->Add(hPMTRateBG,-1.0);
        
        fit_pks(hPMTRate, pmt_res[i], pmt_n0[i]);
        char res_lbl[1024];
        sprintf(res_lbl,"%.0f chn/PE",1000./pmt_res[i]);
        TLatex L(9,hPMTRate->GetMaximum(),res_lbl);
        char light_lbl[1024];
        sprintf(light_lbl,"%.1f PE/MeV",pmt_n0[i]/0.9948);
        TLatex L2(9,0.87*hPMTRate->GetMaximum(),light_lbl);
            
        //hPMTRate->SetMaximum(0.25);
        hPMTRate->Draw();
        hPMTRateBG->Draw("Same");
        L.Draw();
        L2.Draw();
        myA->printCanvas("Signal_"+itos(i+1));
    }
    
    double total_PE = 0;
    printf("PMT response, Channels per PE:\n");
    for(unsigned int i=0; i<N_E_PMT; i++) {
        double sigPerPE = 1000./pmt_res[i];
        double PEperMeV = pmt_n0[i]/0.9948;
        printf("\t%i)\t%g\t%g\n",i, sigPerPE, PEperMeV);
        total_PE += pmt_n0[i]/0.9948;
        
        double sigPerMeV = sigPerPE * PEperMeV;
        AcornDB::ADB().loadPMTcal(RunID(1325,0), RunID(1331,9999), i, 1000./pmt_res[i], sigPerMeV);
    }
    printf("Total %g PE/MeV\n\n", total_PE);
}
