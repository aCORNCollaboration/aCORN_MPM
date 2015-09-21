/// \file ResponseMatrix.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "ResponseMatrix.hh"
#include <cmath>
#include <cassert>
#include <vector>
using std::vector;

void ResponseMatrix::gaussian_response(int inbin, double mu, double sigma, double a) {
    if(inbin < 0 || inbin >= M.GetNrows()) return;
    double n = 1./(sqrt(2*M_PI)*sigma);
    for(int c = 0; c <  M.GetNcols(); c++) {
        double x = out_center(c);
        M(inbin,c) += a*n*exp(-(x-mu)*(x-mu)/(2*sigma*sigma));
    }
}

void ResponseMatrix::normalize_columns() {
    for(int c = 0; c <  M.GetNcols(); c++) {
        double s = 0;
        for(int r = 0; r <  M.GetNrows(); r++) s += M(r,c);
        if(!s) continue;
        for(int r = 0; r <  M.GetNrows(); r++) M(r,c) /= s;
    }
}

void ResponseMatrix::apply(const float* in, float* out) const {
    TMatrix Vin;
    Vin.Use(M.GetNcols(),1,in);
    TMatrix Vout;
    Vout.Use(M.GetNrows(),1,out);
    Vout.Mult(M,Vin);
}

void ResponseMatrix::apply(const TH1F& hIn, TH1F& hOut) const {
    assert(hIn.GetNbinsX() == M.GetNrows());
    assert(hOut.GetNbinsX() == M.GetNcols());
    apply(hIn.GetArray()+1, hOut.GetArray()+1);
}

void ResponseMatrix::smearColumns(const TH2& hIn, TH2& hOut, double s) const {
    assert(hIn.GetNbinsY() == M.GetNrows());
    assert(hOut.GetNbinsY() == M.GetNcols());
    
    vector<float> v_in(M.GetNrows());
    vector<float> v_out(M.GetNcols());
        
    int nx = hIn.GetNbinsX();
    assert(nx == hOut.GetNbinsX());
    for(int ix = 1; ix <= nx; ix++) {
        for(size_t iy=0; iy<v_in.size(); iy++) v_in[iy] = hIn.GetBinContent(hIn.GetBin(ix,iy+1));
        apply(v_in.data(), v_out.data());
        for(size_t iy=0; iy<v_out.size(); iy++) hOut.AddBinContent(hOut.GetBin(ix,iy+1), s*v_out[iy]);

    }
}
