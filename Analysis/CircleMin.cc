/// \file CircleMin.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2016

#include "CircleMin.hh"
#include <cmath>
#include <iostream>
using std::cout;

#include <Math/Functor.h>

double CircleMin::circleMin(const double* params) {
    double s_err = 0;
    
    iSigma(0,0) = params[2];
    iSigma(0,1) = iSigma(1,0) = params[3];
    iSigma(1,1) = params[4];
    iSigma.invert();
    
    for(size_t i=0; i<xs.size(); i++) {
        cs[i] = xs[i]-params[0];
        ss[i] = ys[i]-params[1];
        rs[i] = sqrt(cs[i]*cs[i]+ss[i]*ss[i]);
        cs[i] /= rs[i];
        ss[i] /= rs[i];
        double r2e = (iSigma(0,0)*cs[i]+iSigma(0,1)*ss[i])*cs[i] + (iSigma(1,0)*cs[i]+iSigma(1,1)*ss[i])*ss[i];
        rfits[i] = 1./sqrt(r2e);
        if(!(r2e > 0)) rfits[i] = 0;
        s_err += pow(rs[i]-rfits[i],2);
    }
    //cout << params[0] << "\t" <<  params[1] << "\t" <<  params[2] << "\t" <<  params[4] << "\t" <<  s_err << "\n";
    assert(s_err == s_err);
    return s_err/xs.size();
}

void CircleMin::initGuess(double& x0, double& y0, double& r0) {
    size_t npts = xs.size();
    rs.resize(npts);
    rfits.resize(npts);
    cs.resize(npts);
    ss.resize(npts);
    
    x0 = y0 = 0;
    for(size_t i=0; i<npts; i++) {
        x0 += xs[i];
        y0 += ys[i];
    }
    x0 /= npts;
    y0 /= npts;
    double params[5] = {x0, y0, 1e-6, 0, 1e-6};
    r0 = sqrt(circleMin(params));
}

TGraph* CircleMin::ptsGraph() const {
    TGraph* g = new TGraph(xs.size());
    for(size_t i=0; i<xs.size(); i++)
        g->SetPoint(i,xs[i],ys[i]);
    return g;
}

double CircleMin::doFit() {

    min.SetMaxFunctionCalls(1000);
    min.SetMaxIterations(1000);
    min.SetTolerance(0.0001);
    
    const unsigned int nvar = 5;
    ROOT::Math::Functor f(this, &CircleMin::circleMin, nvar);
    double x0, y0, r0;
    initGuess(x0,y0,r0);
    double variable[nvar] = {x0, y0, r0*r0, 0.001*r0*r0, r0*r0 };
    double step[nvar] = {r0/10., r0/10., variable[2]/10., variable[2]/10., variable[2]/10.};
    cout << "Initial guess: x = " << x0 << " y = " << y0 << " r = " << r0 << "\n";
    min.SetFunction(f);
    
    // Set the free variables to be minimized!
    min.SetVariable(0,"x",variable[0], step[0]);
    min.SetVariable(1,"y",variable[1], step[1]);
    min.SetLimitedVariable(2,"rxx",variable[2], step[2], 0, 1.5*variable[2]);
    min.SetLimitedVariable(3,"rxy",variable[3], step[3], 0, 1.5*variable[2]);
    min.SetLimitedVariable(4,"ryy",variable[4], step[4], 0, 1.5*variable[2]);
    
    min.Minimize();
    
    const double *xs = min.X();
    verbose = true;
    double rms = sqrt(circleMin(xs));
    cout << "Minimum: f( ";
    for(unsigned int i=0; i<nvar; i++) cout << xs[i] << " ";
    cout << "): rms = " << rms << "\n";
    return rms;
}

void CircleMin::transform(double x0, double y0, const Matrix<2,2,double>& M) {
    for(size_t i=0; i<xs.size(); i++) {
        xs[i] -= x0;
        ys[i] -= y0;
        double xnew = M(0,0)*xs[i] + M(0,1)*ys[i];
        ys[i] = M(1,0)*xs[i] + M(1,1)*ys[i];
        xs[i] = xnew;
    }
}

