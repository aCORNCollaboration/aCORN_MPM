#include "HistogramSequenceFitter.hh"
#include <TVirtualFitter.h>
#include <TGraphErrors.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <cassert>
#include <stdio.h>

void HistogramSequenceFitter::addData(const TH1* h, double t1, double t2) {
	assert(h);
	if(hs.size()) assert(h->GetNbinsX() == hs[0]->GetNbinsX()
						 && h->GetNbinsY() == hs[0]->GetNbinsY()
						 && h->GetNbinsZ() == hs[0]->GetNbinsZ());
	
	hs.push_back(h);
	T1.push_back(t1);
	T2.push_back(t2);
}

void HistogramSequenceFitter::fit() {
	assert(hs.size());
	
	TGraph2DErrors g(hs.size());
	fts.clear();
	
	unsigned int nbins = dynamic_cast<const TArray*>(hs[0])->GetSize();
	printf("Fitting sequence of %i histograms with %i bins...\n", (int)hs.size(), nbins);
	for(unsigned int i=1; i<nbins-1; i++) {
		// collect data
		for(unsigned int j=0; j<hs.size(); j++) {
			g.SetPoint(j, T1[j], T2[j], hs[j]->GetBinContent(i));
			g.SetPointError(j, 0, 0, hs[j]->GetBinError(i));
		}
		
		// fit
		estimate_params(g);
		//printf("Parameter estimates:"); for(int j=0; j<f->GetNpar(); j++) printf("\t%g",f->GetParameter(j));
		fts.push_back(g.Fit(f,"S Q N 0"));
		//printf("\nResulting parameter:"); for(int j=0; j<f->GetNpar(); j++) printf("\t%g",f->GetParameter(j));
		//printf("\n");
	}
}

TH1* HistogramSequenceFitter::interpolate(double t1, double t2) {
	assert(hs.size() && fts.size());
	TH1* h = (TH1*)hs[0]->Clone();
	
	double x[2] = {t1, t2};
	unsigned int npar = f->GetNpar();
	TVectorD ig(npar);
	
	for(unsigned int i=0; i<fts.size(); i++) {
		f->SetFitResult(*(fts[i]));
		h->SetBinContent(i+1, (*f)(x));
		// error calculation
		TMatrixDSym covMatrix = fts[i]->GetCovarianceMatrix();
		for(unsigned int j=0; j<npar; j++) {
			if(covMatrix(j,j) > 0 ) {
				ig[j] = f->GradientPar(j,x);
			} else {
				ig[j] = 0; // skip parameters with 0 error
			}
		}
		h->SetBinError(i+1, sqrt(covMatrix.Similarity(ig)));
	}
	
	return h;
}

//----------------

TF1 ExponentialSequenceFitter::expFit = TF1("ExponentialSequenceFitter_Estimator","expo",0,1);

void ExponentialSequenceFitter::estimate_params(const TGraph2DErrors& d) {
	TGraphErrors g(d.GetN());
	for(int i=0; i<d.GetN(); i++) {
		double dt = d.GetY()[i] - d.GetX()[i];
		g.SetPoint(i, (d.GetX()[i] + d.GetY()[i])*0.5 - T0, d.GetZ()[i]/dt);
		g.SetPointError(i, 0, d.GetEZ()[i]/dt);
	}
	g.Fit(&expFit,"QN0");
	f->SetParameter(0,exp(expFit.GetParameter(0))/expFit.GetParameter(1));
	f->SetParameter(1,expFit.GetParameter(1));
}

//----------------

double PolynomialSequenceFitter::integ_f(double t, const double* params) const {
	double x = t-T0;
	double xn = 1;
	double s = 0;
	for(int i=0; i<f->GetNpar(); i++) {
		xn *= x;
		s += params[i]*xn/(i+1);
	}
	return s;
}

void PolynomialSequenceFitter::estimate_params(const TGraph2DErrors& d) {
	TGraphErrors g(d.GetN());
	for(int i=0; i<d.GetN(); i++) {
		double dt = d.GetY()[i] - d.GetX()[i];
		g.SetPoint(i, (d.GetX()[i] + d.GetY()[i])*0.5 - T0, d.GetZ()[i]/dt);
		g.SetPointError(i, 0, d.GetEZ()[i]/dt);
	}
	char c[128];
	sprintf(c,"pol%i",f->GetNpar()-1);
	TF1 polFit("polFit",c,0,1);
	g.Fit(&polFit,"QN0");
	for(int i=0; i<f->GetNpar(); i++) f->SetParameter(i, polFit.GetParameter(i));
}
