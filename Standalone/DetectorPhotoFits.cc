/// \file DetectorPhotoFits.cc Fit circles on photo of masked detector
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

// make DetectorPhotoFits
// ./DetectorPhotoFits /home/mpmendenhall/Downloads/Proton\ Mask\ Photos/P1030627.JPG_points.txt

#include "CircleMin.hh"
#include "StringManip.hh"
#include "GraphicsUtils.hh"
#include "PathUtils.hh"

#include <TSystem.h>
#include <TPad.h>
#include <TAxis.h>
#include <cmath>
#include <fstream>
#include <streambuf>
#include <TMarker.h>
#include <TLine.h>
#include <TLatex.h>

#include <cmath>
#include <iostream>
using std::cout;


void DetectorCircles(const string& datfname) {
    
    vector<CircleMin*> circles;
    
    vector<string> datlines = split(loadFileString(datfname),"\n");
    int gprev = -1;
    for(auto it = datlines.begin(); it != datlines.end(); it++) {
        vector<double> p = sToDoubles(*it,"\t");
        if(p.size() != 3) break;
        if(p[0] != gprev) {
            circles.push_back(new CircleMin());
            gprev = p[0];
        }
        circles.back()->addPoint(p[1],p[2]);
    }
    
    printf("Loaded %zu circles.\n", circles.size());
    
    int linecolor = 2;
    
    circles[0]->doFit();
    const double* xs0 = circles[0]->min.X();
    double x0 = xs0[0];
    double y0 = xs0[1];
    auto M = circles[0]->iSigma;
    double r0 = 4.5/2.;
    M(0,0) = r0*sqrt(M(0,0));
    M(1,1) = r0*sqrt(M(1,1));
    M(0,1) = M(1,0) = 0;
    
    for(auto C: circles) {
        C->transform(x0,y0,M);
        
        C->doFit();
        const double *xs = C->min.X();
        
        TGraph* g = C->ptsGraph();
        g->SetMarkerStyle(7);
        g->Draw(C==circles[0]? "AP" : "P");
        g->SetTitle("image circle fits");
        g->GetXaxis()->SetTitle("horizontal [a.u.]");
        g->GetYaxis()->SetTitle("vertical [a.u.]");
        g->GetYaxis()->SetTitleOffset(1.5);
        g->GetXaxis()->SetLimits(-3,3);
        g->GetYaxis()->SetRangeUser(-3,3);
        
        TPolyLine* LE = makeEllipse(xs[0], xs[1], &C->iSigma[0]);
        LE->SetLineColor(linecolor);
        LE->Draw();
        
        if(C != circles[0]) {
            char lbl[1024];
            //sprintf(lbl,"%+.1f, %+.1f", 1000*xs[0], 1000*xs[1]); // centers, in mils
            double r_m = pow(xs[2]*xs[4],0.25);
            sprintf(lbl,"%.3f", r_m); // radii
            
            TLatex* L = new TLatex(xs[0]+0.05,xs[1]+0.05,lbl);
            L->SetTextColor(linecolor);
            L->SetTextSize(0.025);
            L->Draw();
            
            printf("r1 = %g, r2 = %g, xy = %g\n", sqrt(xs[2]), sqrt(xs[4]), sqrt(xs[3]));
            double d = 4.5;
            double r_b = 1.725;
            double r_t = 2.25;
            double u = 1 - r_b/r_t;            
            double r_s = d*(u-1 + sqrt((1-u)*(1-u) + 4*u*r_m/d))/(2*u);
            printf("u = %g, r_m = %g: r_s = %g\n", u, r_m, r_s);
        }
        
        // center markers
        TMarker* mCenter = new TMarker(xs[0],xs[1],2);
        mCenter->SetMarkerColor(linecolor);
        mCenter->Draw();
        
        linecolor += 2;
    }
    
    gPad->SetCanvasSize(300,300);
    string outname = dropLast(datfname,".")+".pdf";
    gPad->Print(outname.c_str());
}

int main(int argc, char** argv) {
    if(argc != 2) return EXIT_FAILURE;
    DetectorCircles(argv[1]);
    return EXIT_SUCCESS;
}
