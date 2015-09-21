/// \file SimWishboneSmearer.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "SimWishboneSmearer.hh"
#include <TAxis.h>
#include <cassert>

SimWishboneSmearer::SimWishboneSmearer(const TH2& h, int ybinsOut):
ResponseMatrix(ybinsOut < 0? h.GetNbinsY() : ybinsOut, h.GetNbinsY()), hIn(h) {
    const TAxis* Ay = hIn.GetYaxis();
    in_range[0] = Ay->GetBinLowEdge(1);
    in_range[1] = Ay->GetBinUpEdge(Ay->GetNbins());
}

void SimWishboneSmearer::setResponse(const TH2& hOut) {
    const TAxis* Ay = hOut.GetYaxis();
    assert(Ay->GetNbins() == M.GetNrows());
    out_range[0] = Ay->GetBinLowEdge(1);
    out_range[1] = Ay->GetBinUpEdge(Ay->GetNbins());
    
    for(int iy=0; iy<hIn.GetNbinsY(); iy++) {
        double t = in_center(iy);
        gaussian_response(iy, t*tscale + t0, sigma);
    }
    normalize_columns();
}
