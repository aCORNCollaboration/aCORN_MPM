#include "HistogramSequenceFitter.hh"
#include <TVirtualFitter.h>
#include <TGraphErrors.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <cassert>
#include <stdio.h>


IntervalIntegralFitter::IntervalIntegralFitter(unsigned int npar):
N(npar), 
//myMin(ROOT::Math::kConjugateFR),
myMin(ROOT::Math::kVectorBFGS2),
myf(this, &IntervalIntegralFitter::eval_error, npar) {
    myParams = new double[N];
    mySteps = new double[N];
    myMin.SetFunction(myf);
}

IntervalIntegralFitter::~IntervalIntegralFitter() {
    delete[] myParams;
    delete[] mySteps;
}

void IntervalIntegralFitter::fit() {
    myMin.SetMaxFunctionCalls(10000);
    myMin.SetMaxIterations(10000);
    myMin.SetTolerance(0.1);
    
    const char* varletters = "xyztuvwabcdefghijklmnopqrs";
    for(unsigned int i=0; i<N; i++) {
        const char vname[] = {varletters[i],0};
        //printf("%s: %g\t%g\n", vname, myParams[i], mySteps[i]);
        myMin.SetVariable(i, vname, myParams[i], mySteps[i]);
    }
    myMin.SetVariableLimits(1,-1e-3,1e-3);
    
    if(dIntegrals.size() != integrals.size()) {
        dIntegrals = vector<double>(integrals.size(), 1.0);
    }
    
    myMin.Minimize(); 
    
    const double* xs = myMin.X();
    for(unsigned int i=0; i<N; i++) {
        myParams[i] = xs[i];
        //printf("\t\t=> p%i = %g\n",i,myParams[i]);
    }
}

double IntervalIntegralFitter::eval_error(const double* params) const {
    assert(intervals.size() == integrals.size());
    assert(intervals.size() == dIntegrals.size());
    //printf("p ="); for(unsigned int i=0; i<N; i++) printf("\t%g",params[i]); printf("\n");
    
    double ss = 0;
    for(size_t i = 0; i<intervals.size(); i++) {
        if(!dIntegrals[i]) continue;
        double yhat = (*this)(intervals[i],params);
        double dy = (yhat-integrals[i])/dIntegrals[i];
        //printf("%i) %g(%g)\t%g\n", i, integrals[i], dIntegrals[i], yhat);
        ss += dy*dy;
    }
    return ss;
}

double IntervalIntegralFitter::operator()(const intervalList& L, const double* params) const {
    if(!params) params = myParams;
    double I = 0;
    for(auto it = L.begin(); it != L.end(); it++) I += integ_f(it->second,params) - integ_f(it->first,params);
    return I;
}

//----------------

double mean_point(const intervalList& L, double& w) {
    double sx = 0;
    w = 0;
    for(auto it = L.begin(); it != L.end(); it++) {
        sx += 0.5*(it->second+it->first)*(it->second-it->first);
        w += it->second - it->first;
    }
    return sx/w;
}

//----------------

TF1 ExponentialIntegralFitter::expFit = TF1("ExponentialSequenceFitter_Estimator","expo",0,10000);

double ExponentialIntegralFitter::integ_f(double t, const double* params) const {
    double x = t-T0;
    if(fabs(x*params[1]) < 1e-3) {
        return params[0]*x*(1 + params[1]*x/2. + pow(params[1]*x,2)/6. + pow(params[1]*x,3)/24);
    } else {
        return params[0]/params[1] * exp(x*params[1]) - params[0]/params[1];
    }
}

void ExponentialIntegralFitter::init_params() {    
    //printf("Estimating initial exponential parameters...\n");
    TGraphErrors g(intervals.size());
    for(size_t i=0; i<intervals.size(); i++) {
        double w;
        double x0 = mean_point(intervals[i],w);
        g.SetPoint(i, x0-T0, integrals[i]/w);
        g.SetPointError(i, 0, dIntegrals[i]/w);
        //printf("\t%i] %g (%g)\t%g (%g)\n", i, x0-T0, w, integrals[i]/w, dIntegrals[i]/w);
        myParams[0] = integrals[i]/w;
        myParams[1] = 0;
    }
    //g.Fit(&expFit,"MN");
    //for(int i=0; i<2; i++) {
    //    myParams[i] = expFit.GetParameter(i);
    //    printf("\tp%i = %g\n",i,myParams[i]);
    //}
    for(int i=0; i<2; i++) mySteps[i] = fabs(myParams[i]*0.1);
    mySteps[1] = 1e-3;
}

//----------------

double PolynomialIntegralFitter::integ_f(double t, const double* params) const {
    double x = t-T0;
    double xn = 1;
    double s = 0;
    for(unsigned int i=0; i<N; i++) {
        xn *= x;
        s += params[i]*xn/(i+1);
    }
    return s;
}

void PolynomialIntegralFitter::init_params() {
    /*TGraphErrors g(d.GetN());
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
    */
}




///////////////////////////////////
///////////////////////////////////
///////////////////////////////////


void HistogramSequenceFitter::addData(const TH1* h, const intervalList& dt) {
    assert(h);
    if(hs.size()) assert(h->GetNbinsX() == hs[0]->GetNbinsX()
                         && h->GetNbinsY() == hs[0]->GetNbinsY()
                         && h->GetNbinsZ() == hs[0]->GetNbinsZ());
    hs.push_back(h);
    dts.push_back(dt);
}

void HistogramSequenceFitter::fit() {
    assert(myFitter);
    assert(hs.size() >= myFitter->N);
    
    myFitter->intervals = dts;
    myFitter->integrals.resize(hs.size());
    myFitter->dIntegrals.resize(hs.size());
    fts.clear();
    
    unsigned int nbins = dynamic_cast<const TArray*>(hs[0])->GetSize();
    printf("Fitting sequence of %i histograms with %i bins...\n", (int)hs.size(), nbins);
    for(unsigned int i=0; i<=nbins; i++) {
        // collect data
        for(unsigned int j=0; j<hs.size(); j++) {
            myFitter->integrals[j] = hs[j]->GetBinContent(i);
            myFitter->dIntegrals[j] = hs[j]->GetBinError(i);
        }
        myFitter->init_params();
        myFitter->fit();
        fts.push_back(vector<double>(myFitter->myParams, myFitter->myParams + myFitter->N));
    }
}

TH1* HistogramSequenceFitter::interpolate(const intervalList& dt) {
    assert(hs.size() && fts.size());
    TH1* h = (TH1*)hs[0]->Clone();
    
    //TVectorD ig(npar);
    
    for(unsigned int i=0; i<fts.size(); i++) {
        h->SetBinContent(i, (*myFitter)(dt, &fts[i][0]));
        
        // TODO error calculation
        //TMatrixDSym covMatrix = fts[i]->GetCovarianceMatrix();
        //for(unsigned int j=0; j<npar; j++) {
        //    if(covMatrix(j,j) > 0 ) {
        //        ig[j] = f->GradientPar(j,x);
        //    } else {
        //        ig[j] = 0; // skip parameters with 0 error
        //    }
        //}
        //h->SetBinError(i+1, sqrt(covMatrix.Similarity(ig)));
    
    }
    
    return h;
}

