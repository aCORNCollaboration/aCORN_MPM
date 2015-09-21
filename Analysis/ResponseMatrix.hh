/// \file ResponseMatrix.hh Histogram-smearing response response matrix class
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef RESPONSEMATRIX_HH
#define RESPONSEMATRIX_HH

#include <TMatrix.h>
#include <TH1.h>
#include <TH2.h>

class ResponseMatrix {
public:
    /// Constructor, with rows (output) and columns (input)
    ResponseMatrix(int nr, int nc): M(nr, nc) { }
    
    double in_range[2] = {0,1};         ///< input axis range
    /// center of numbered input bin
    double in_center(int b) const { return in_range[0] + (b+0.5)*(in_range[1]-in_range[0])/M.GetNrows(); }
    
    double out_range[2] = {0,1};        ///< output axis range
    /// center of numbered input bin
    double out_center(int b) const { return out_range[0] + (b+0.5)*(out_range[1]-out_range[0])/M.GetNcols(); }
    
    /// map specified input bin to gaussian on output (area approximately s up to binning normalization error)
    void gaussian_response(int inbin, double mu, double sigma, double a = 1.);
    /// Normalize columns to unit sum
    void normalize_columns();
    
    /// apply to float data
    void apply(const float* in, float* out) const;
    /// apply between TH1Fs
    void apply(const TH1F& hIn, TH1F& hOut) const;
    /// apply to TH2 columns, summing in to output with scale s
    void smearColumns(const TH2& hIn, TH2& hOut, double s = 1.) const;
        
protected:
    TMatrix M;  ///< the matrix
};

#endif
