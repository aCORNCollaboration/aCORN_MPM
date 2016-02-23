/// \file WishboneTemplateRate.cc Low-statistics wishbone rate extraction via template inner product
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

// make WishboneTemplateRate -j4
// ./WishboneTemplateRate

#include "PathUtils.hh"
#include "StringManip.hh"
#include "WishbonePlugin.hh"
#include <TH2.h>


double TH2_dot(const TH2& h0, const TH2& h1, double x0, double x1, double y0, double y1, double& serr) {
    Int_t nc = h0.GetNcells();
    Int_t bx,by,bz;
    double s = 0;
    serr = 0;
    for(int c=0; c<nc; c++) {
        h1.GetBinXYZ(c,bx,by,bz);
        double x = h1.GetXaxis()->GetBinCenter(bx);
        double y = h1.GetYaxis()->GetBinCenter(by);
        if(!(x0 <= x && x < x1 && y0 <= y && y < y1)) continue;
        
        double dz = h0.GetBinError(c);
        double dz2 = h1.GetBinError(c);
        double w = fabs(h0.GetBinContent(c)/pow(dz,2)); // near-optimal
        //double w = 1; // awful
            
        s += w * h1.GetBinContent(c);
        serr += pow(w * dz2, 2);
    }
    serr = sqrt(serr);
    printf("%g +- %g    (%g)\n", s, serr, serr/s);
    return s;
}

int main(int, char**) {
    
    OutputManager OM("NameUnused", "/home/mpmendenhall/Desktop/foo");
    string bpath = getEnvSafe("ACORN_WISHBONE")+"/Series_";
    WishboneAnalyzer WA0(&OM, "NameUnused", bpath+"3532/Series_3532.root");
    WishboneAnalyzer WA1(&OM, "NameUnused", bpath+"3540/Series_3540.root");
    WishboneAnalyzer WA2(&OM, "NameUnused", bpath+"3541/Series_3541.root");
    WA0.calculateResults();
    WA1.calculateResults();
    WA2.calculateResults();
    
    TH2& h0 = *((WishbonePlugin*)WA0.getPlugin("WishbonePlugin"))->hWishboneBGSub;
    TH2& h1 = *((WishbonePlugin*)WA1.getPlugin("WishbonePlugin"))->hWishboneBGSub;
    TH2& h2 = *((WishbonePlugin*)WA2.getPlugin("WishbonePlugin"))->hWishboneBGSub;
    double ds0, ds1, ds2;
    double t0 = 3.0, t1 = 4.5;
    double d0 = TH2_dot(h0,h0, 50, 600, t0, t1, ds0);
    double d1 = TH2_dot(h0,h1, 50, 600, t0, t1, ds1);
    double d2 = TH2_dot(h0,h2, 50, 600, t0, t1, ds2);
    printf("d0 = %g; d1 = %g; d1/d0 = %g (%g);  d2/d0 = %g (%g)\n", d0, d1, d1/d0, ds1/d0, d2/d0, ds2/d0);
    
    AcornDB::closeDB();
    return EXIT_SUCCESS;
}
