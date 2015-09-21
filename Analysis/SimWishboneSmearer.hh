/// \file SimWishboneSmearer.hh Routine for smearing/scaling simulated wishbone for global fit
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef SIMWISHBONESMEARER_HH
#define SIMWISHBONESMEARER_HH

#include "ResponseMatrix.hh"

class SimWishboneSmearer: protected ResponseMatrix {
public:
    /// Constructor
    SimWishboneSmearer(const TH2& h, int ybinsOut = -1);
    
    TH2 const& hIn;     ///< input simulated wishbone
    
    /// set up response matrix for output histogram binning
    void setResponse(const TH2& hOut);
    /// sum into output histogram
    void sumOutput(TH2& hOut) { smearColumns(hIn, hOut, zscale); }
    
    double sigma = 1.0; ///< Gaussian time smearing sigma [us]
    double t0 = 0;      ///< time axis offset; t -> t*tscale + t0
    double tscale = 1;  ///< time axis scaling; t -> t*tscale + t0
    double zscale = 1;  ///< z scaling
};

#endif
