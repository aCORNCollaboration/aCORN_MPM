/// \file CircleFit.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

// make CircleFit; ./CircleFit; display circle.pdf

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

string m2dms(double m, bool texmode = false) {
    int dg = floor(m/60);
    m -= 60*dg;
    int mn = floor(m);
    double sc = 60*(m-mn);
    if(texmode) return to_str(dg) + "#circ" + to_str(mn) + "'" + to_str(int(sc)) + "\"";
    return to_str(dg) + "° " + to_str(mn) + "' " + to_str(int(sc)) + "\"";
}

void do_circle_fit(const string& fbase) {
    CircleMin QM;
    std::ifstream t((fbase+".csv").c_str());
    string dat((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());
    cout << "------------ Input file:\n" << dat << "-----------------\n\n";
    vector<string> datlines = split(dat,"\n");
    double cxhx[2] = {0,0};
    for(auto it = datlines.begin(); it != datlines.end(); it++) {
        vector<double> pts = sToDoubles(*it,",");
        vector<string> wds = split(*it,",");
        if(pts.size() >= 12) {
            if(wds[9] == "Crosshair") {
                if(pts.size()<18) for(int i=0; i<2; i++) cxhx[i] = pts[10+i]/60.;
                else for(int i=0; i<2; i++) cxhx[i] = pts[16+i]/60.;
                cxhx[1] = 90*60-cxhx[1];
            }
        }
        if(pts.size() < 10 || !pts[7]) continue;
        cout << pts.size() << ") '" << strip(*it) << "':\t" << pts[7] << "\t" << pts[8] << "\n";
        QM.addPoint(pts[7]/60.,90*60-pts[8]/60.);
    }
    
    double rms = QM.doFit();
    const double *xs = QM.min.X();

    TGraph* g = QM.ptsGraph();
    g->SetTitle("theodolite circle fit");
    g->GetXaxis()->SetTitle("horizontal [arcmin]");
    g->GetYaxis()->SetTitle("90#circ-vertical [arcmin]");
    g->GetYaxis()->SetTitleOffset(1.5);
    g->Draw("A*");
    
    TGraph gCirc(100);
    double rr = sqrt(0.5*(xs[2]+xs[4]));
    if(fabs(rr-sqrt(xs[2])) < 0.05*rr) {
        for(int i=0; i<100; i++) {
            double th = i*2*M_PI/99.;
            gCirc.SetPoint(i, xs[0]+rr*cos(th), xs[1]+rr*sin(th));
        }
        gCirc.SetLineColor(4);
        gCirc.Draw("C");
    }
        
    TPolyLine* LE = makeEllipse(xs[0], xs[1], &QM.iSigma[0]);
    LE->SetLineColor(2);
    LE->Draw();
    
    // exaggerated errors
    if(rms < 0.005*rr) {
        for(size_t i=0; i<QM.xs.size(); i++) {
            double r = QM.rfits[i]+100*(QM.rs[i]-QM.rfits[i]);
            double x1 = xs[0] + QM.cs[i]*r;
            double y1 = xs[1] + QM.ss[i]*r;
            TMarker* mOff = new TMarker(x1, y1, 4);
            mOff->SetMarkerColor(4);
            mOff->SetMarkerSize(0.5);
            mOff->Draw();
            TLine* lOff = new TLine(QM.xs[i], QM.ys[i], x1, y1);
            lOff->SetLineColor(4);
            lOff->Draw();
        }
    }
    
    // center markers
    TMarker mCenter(xs[0],xs[1],2);
    mCenter.SetMarkerColor(2);
    mCenter.Draw();
    printf("Center at\tH = %s\tV = %s\n",m2dms(xs[0]).c_str(),m2dms(90*60-xs[1]).c_str());

    char lbl[1024];
    sprintf(lbl,"#splitline{Center (%s, %s), RMS %.2f'}{r_{xx} = %.2f'; r_{yy} = %.2f'; r_{xy} = %.2f}",
            m2dms(xs[0],true).c_str(),m2dms(90*60-xs[1],true).c_str(), rms, sqrt(xs[2]), sqrt(xs[4]), sqrt(fabs(xs[3])));
    TLatex L(xs[0]-0.7*rr,xs[1]+0.25*rr,lbl);
    L.SetTextSize(0.03);
    L.Draw();

    double cxdist = sqrt(pow(xs[0]-cxhx[0],2)+pow(xs[1]-cxhx[1],2));
    if(cxdist < sqrt(xs[2])) {
        TMarker* mCx = new TMarker(cxhx[0],cxhx[1],5);
        mCx->SetMarkerColor(4);
        mCx->Draw();
        sprintf(lbl,"Crosshair offset: (%.2f', %.2f')", cxhx[0]-xs[0], xs[1]-cxhx[1]);
        TLatex* L2 = new TLatex(xs[0]-0.6*rr,xs[1]-0.4*rr,lbl);
        L2->SetTextSize(0.03);
        L2->Draw();
        printf("Crosshair at\tH = %s\tV = %s\n",m2dms(cxhx[0]).c_str(),m2dms(90*60-cxhx[1]).c_str());
    }
    
    gPad->SetCanvasSize(300,300);
    gPad->Print((fbase+".pdf").c_str());
}

int main(int argc, char** argv) {
    string basedir = (argc>1)?argv[1]:"/home/mpmendenhall/Documents/aCORN/20150112_Alignment";
    vector<string> fnames = listdir(basedir);
    for(auto it = fnames.begin(); it != fnames.end(); it++) {
        if(split(*it," ")[0] != "Finding" || split(*it,".").back() != "csv") continue;
        do_circle_fit(basedir+"/"+dropLast(*it,"."));
    }
}
