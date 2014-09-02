#include "ReducedDataScanner.hh"
#include "OutputManager.hh"
#include "strutils.hh"
#include "HistogramSequenceFitter.hh"
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>

double fill_energy_histogram(ReducedDataScanner& R, TH1* h, double escale = 1000./38000.) {
	R.startScan();
	h->Sumw2();
	Long_t clk_ttl = 0;
	Long_t clk_prev = 0;
	while(R.nextPoint()) {
		if(R.T_p < clk_prev) clk_ttl += clk_prev;
		clk_prev = R.T_p;
	
		if(R.nV) continue;
		h->Fill(R.E_e * escale);
	}
	h->GetXaxis()->SetTitle("energy [keV]");
	//h->SetMinimum(1);
	//h->SetMaximum(1e6);
	return (clk_prev+clk_ttl)*(1e-8);
}

void bg_subtraction_study() {

	OutputManager OM("Bi_PMT_Gain","/Users/michael/Documents/aCORN_plots/");
	gStyle->SetOptStat("");
	OM.defaultCanvas->SetLogy(true);
	OM.defaultCanvas->SetLeftMargin(0.14);
	
	double nbins = 200;
	double emax = 200;
	double yscale = 1.0;
	const char* ylbl = "event rate [Hz/keV]";
	
	std::vector<ReducedDataScanner*> Rb;
	std::vector<TH1F*> hBG;
	//ExponentialSequenceFitter ESF;
	PolynomialSequenceFitter ESF(2);
	
	for(unsigned int i=0; i<3; i++) {
		Rb.push_back(new ReducedDataScanner);
		Rb.back()->addFile("/Users/michael/Documents/aCORN_data/ROOTified/s"+itos(1326+2*i)+"r00*");
		hBG.push_back(new TH1F(("hBG_"+itos(i)).c_str(),"Source events spectrum",nbins,0,emax));
		double dt = fill_energy_histogram(*Rb.back(), hBG.back());
		printf("Loaded %.1fs background.\n",dt);
		hBG.back()->Draw();
		OM.printCanvas("BG_spectrum_"+itos(i));
		
		ESF.addData(hBG.back(), i*60*60, i*60*60 + dt);
	}
	
	ESF.fit();
	

	
	ReducedDataScanner Rf;
	Rf.addFile("/Users/michael/Documents/aCORN_data/ROOTified/s1327r00*");
	TH1F* hBi = new TH1F("hBi","Bi207 events spectrum",nbins,0,emax);
	double dt = fill_energy_histogram(Rf,hBi);
	printf("Loaded %.1fs foreground.\n",dt);
	TH1* hBiBg = ESF.interpolate(20*60, 20*60+dt);
	
	hBi->SetLineColor(2);
	hBi->Draw();
	hBiBg->Draw("Same");
	OM.printCanvas("BG_spectrum_Bi");
	
	hBi->Add(hBiBg,-1.0);
	hBi->Scale(yscale/dt/hBi->GetBinWidth(1));
	hBi->GetYaxis()->SetTitle(ylbl);
	hBi->GetYaxis()->SetTitleOffset(1.5);
	
	hBiBg->Scale(yscale/dt/hBiBg->GetBinWidth(1));
	hBiBg->GetYaxis()->SetTitle(ylbl);
	hBiBg->GetYaxis()->SetTitleOffset(1.5);
	
	OM.defaultCanvas->SetLogy(false);
	
	//hBi->SetMinimum(-100);
	//hBiBg->SetMaximum(80);
	hBi->Draw();
	hBiBg->Draw("Same");
	OM.printCanvas("Bi_spectrum");
}

void collect_hdata(TH1F** hBi, ReducedDataScanner& R, double* p_scale, double* p_weight, double weight_sum, double E_min, double E_max) {
	for(unsigned int i=0; i<N_E_PMT; i++) hBi[i]->Scale(0);

	R.startScan();
	while(R.nextPoint()) {
		if(R.nV) continue;
		double E = 0;
		if(p_scale && p_weight) {
			for(unsigned int i=0; i<N_E_PMT; i++) E += R.E_PMT[N_V_PMT+i]*p_scale[i]*p_weight[i];
			E /= weight_sum;
		} else {
			E = R.E_e / weight_sum;
		}
				
		if(E_min < E && E < E_max) {
			for(unsigned int i=0; i<N_E_PMT; i++)
				if(R.E_PMT[N_V_PMT+i]>10)
					hBi[i]->Fill(R.E_PMT[N_V_PMT+i]*p_scale[i]);
		}
	}
}

Double_t poissonf(Double_t* x, Double_t* par) {
  return par[0]*TMath::Poisson(x[0]*par[1], par[2]);
}

double fit_pks(TH1F** hBi, double* p_scale, double* p_weight, bool show_fit = false) {
	TF1 fGaus("fGaus","gaus",0,100);
	TF1 fPois("pois", &poissonf, 0, 10, 3);
	
	double weight_sum = 0;
	for(unsigned int i=0; i<N_E_PMT; i++) {
		double c0 = hBi[i]->GetBinCenter(hBi[i]->GetMaximumBin());
		fGaus.SetRange(c0-0.3*c0,c0+0.3*c0);
		hBi[i]->Fit(&fGaus,"QRN");
		
		c0 = fGaus.GetParameter(1);
		double w = fGaus.GetParameter(2);
		double N = c0*c0/(w*w);
		
		std::cout << "Estimating N = " << N << " from width w = " << w << " out of c = " << c0 << "\n";
		fPois.SetParameter(0, 1.0);
		fPois.SetParameter(1, N/c0);
		fPois.SetParameter(2, N);
		std::cout << "Initial height = " << fPois.Eval(c0) << "\n";
		fPois.SetParameter(0,fGaus.GetParameter(0)/fPois.Eval(c0));
		fPois.SetRange(c0-2*w,c0+2*w);
		
		hBi[i]->Fit(&fPois,show_fit?"R":"QRN");
		
		p_scale[i] *= 1000./(fPois.GetParameter(2)/fPois.GetParameter(1));
		p_weight[i] = fPois.GetParameter(2);
		weight_sum += p_weight[i];
	}
	return weight_sum;
}

int main(int, char**) {
	
	bg_subtraction_study();
	return 0;

	OutputManager OM("Bi_PMT_Gain","/Users/michael/Documents/aCORN_plots/");
	OM.defaultCanvas->SetLogy(true);
	
	ReducedDataScanner R;
	R.addFile("/Users/michael/Documents/aCORN_data/ROOTified/s1327r00*");

	TH1F* hBi[N_E_PMT];
	TH1F* hBi2[N_E_PMT];
	double p_scale[N_E_PMT];
	double p_weight[N_E_PMT];
	
	for(unsigned int i=0; i<N_E_PMT; i++) {
		hBi[i] = new TH1F(("hBi_"+itos(i)).c_str(),"PMT 207Bi spectrum",100,0,3000);
		hBi2[i] = new TH1F(("hBi2_"+itos(i)).c_str(),"PMT 207Bi spectrum",100,0,3000);
		p_scale[i] = 0.3;
	}
	
	collect_hdata(hBi, R, p_scale, NULL, 38000./1000., 800, 1250);

	// fit and plot
	double weight_sum = 1;
	weight_sum = fit_pks(hBi, p_scale, p_weight, false);
	collect_hdata(hBi,  R, p_scale, NULL, 38000./1000., 800, 1250);
	weight_sum = fit_pks(hBi, p_scale, p_weight, true);
	collect_hdata(hBi2, R, p_scale, NULL, 38000./1000., 0, 3000);
	for(unsigned int i=0; i<N_E_PMT; i++) {
		hBi[i]->SetLineColor(2);
		hBi2[i]->SetLineColor(4);
		hBi2[i]->GetXaxis()->SetTitle("energy [keV]");
		hBi2[i]->GetYaxis()->SetTitle("count");
		hBi2[i]->SetMinimum(2);
		hBi2[i]->Draw();
		hBi[i]->Draw("Same");
		OM.printCanvas("hBi_"+itos(i));
	}
	
	for(unsigned int i=0; i<N_E_PMT; i++) {
		printf("PMT %i:\tscale = %f\tweight = %f\tmul = %f\n",i,p_scale[i],p_weight[i],p_scale[i]*p_weight[i]);
	}
	
	
	weight_sum *= 1.05; // TODO shouldn't need this!
	
	
	
	// re-weighted data
	TH1F* hBiRecon = new TH1F("hBiRecon","PMT 207Bi spectrum",200,0,3000);
	TH1F* hBiScaled = new TH1F("hBiScaled","PMT 207Bi spectrum",200,0,3000);
	R.startScan();
	while(R.nextPoint()) {
		if(R.nV) continue;
		hBiScaled->Fill(R.E_e * 1000./38000.);
		double E = 0;
		for(unsigned int i=0; i<N_E_PMT; i++) E += R.E_PMT[N_V_PMT+i]*p_scale[i]*p_weight[i];
		E /= weight_sum;
		hBiRecon->Fill(E);
	}
	
	hBiRecon->SetLineColor(2);
	hBiScaled->SetLineColor(4);
	
	hBiScaled->Draw();
	hBiRecon->Draw("Same");
	OM.printCanvas("hBi_Reconstructed");
	
	return 0;
}
