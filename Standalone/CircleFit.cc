// make CircleFit; ./CircleFit; display circle.pdf
 
#include "StringManip.hh"
#include "GraphicsUtils.hh"
#include "Matrix.hh"
#include "PathUtils.hh"

#include <Math/GSLMinimizer.h>
#include <Math/Functor.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TPad.h>
#include <TAxis.h>
#include <cmath>
#include <fstream>
#include <streambuf>
#include <TMarker.h>
#include <TLine.h>
#include <TLatex.h>

#include <iostream>
using std::cout;

class CircleMin {
public:
    CircleMin() { }
    
    /// params: x,y,rxx,rxy,ryy
    double circleMin(const double* params) {
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
            rfits[i] = 1./sqrt((iSigma(0,0)*cs[i]+iSigma(0,1)*ss[i])*cs[i] + (iSigma(1,0)*cs[i]+iSigma(1,1)*ss[i])*ss[i]);
            s_err += pow(rs[i]-rfits[i],2);
        }
        return s_err/xs.size();
    }
    
    bool verbose = false;
    vector<double> xs;          ///< x coordinate each point
    vector<double> ys;          ///< y coordinate each point
    vector<double> cs;          ///< cosines each point
    vector<double> ss;          ///< sines each point
    vector<double> rs;          ///< radius each point
    vector<double> rfits;       ///< fit radius each point
    Matrix<2,2,double> iSigma;  ///< inverse covariance matrix
    
    /// calculate initial guess
    void initGuess(double& x0, double& y0, double& r0) {
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
    
    TGraph* ptsGraph() const {
        TGraph* g = new TGraph(xs.size());
        for(size_t i=0; i<xs.size(); i++)
            g->SetPoint(i,xs[i],ys[i]);
        return g;
    }
};

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
            if(wds[9] == "Crosshair")
                for(int i=0; i<2; i++) cxhx[i] = pts[10+i]/60.;
        }
        if(pts.size() < 10 || !pts[7]) continue;
        cout << pts.size() << ") '" << strip(*it) << "':\t" << pts[7] << "\t" << pts[8] << "\n";
        QM.xs.push_back(pts[7]/60.);
        QM.ys.push_back(pts[8]/60.);
    }
    
    // Choose method upon creation between:
    // kConjugateFR, kConjugatePR, kVectorBFGS,
    // kVectorBFGS2, kSteepestDescent
    ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS );
    
    min.SetMaxFunctionCalls(1000);
    min.SetMaxIterations(1000);
    min.SetTolerance(0.0001);
    
    const unsigned int nvar = 5;
    ROOT::Math::Functor f(&QM,&CircleMin::circleMin, nvar);
    double x0, y0, r0;
    QM.initGuess(x0,y0,r0);
    double variable[nvar] = {x0, y0, r0*r0, 0.01*r0*r0, r0*r0 };
    double step[nvar] = {r0/100., r0/100., variable[2]/100., variable[2]/100., variable[2]/100.};
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
    QM.verbose = true;
    double rms = sqrt(QM.circleMin(xs));
    cout << "Minimum: f( ";
    for(unsigned int i=0; i<nvar; i++) cout << xs[i] << " ";
    cout << "): rms = " << rms << "\n";
    
    
    TGraph* g = QM.ptsGraph();
    g->SetTitle("theodolite circle fit");
    g->GetXaxis()->SetTitle("horizontal [arcmin]");
    g->GetYaxis()->SetTitle("vertical [arcmin]");
    g->GetYaxis()->SetTitleOffset(1.5);
    g->Draw("A*");
    
    TGraph gCirc(100);
    double rr = sqrt(0.5*(xs[2]+xs[4]));
    if(abs(xs[2]-xs[4]) < 0.1*rr) {
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
    printf("Center at\tH = %s\tV = %s\n",m2dms(xs[0]).c_str(),m2dms(xs[1]).c_str());

    char lbl[1024];
    sprintf(lbl,"#splitline{Center (%s, %s), RMS %.2f'}{r_{xx} = %.2f'; r_{yy} = %.2f'; r_{xy} = %.2f}",
            m2dms(xs[0],true).c_str(),m2dms(xs[1],true).c_str(), rms, sqrt(xs[2]), sqrt(xs[4]), sqrt(fabs(xs[3])));
    TLatex L(xs[0]-0.7*rr,xs[1]+0.25*rr,lbl);
    L.SetTextSize(0.03);
    L.Draw();

    double cxdist = sqrt(pow(xs[0]-cxhx[0],2)+pow(xs[1]-cxhx[1],2));
    if(cxdist < sqrt(xs[2])) {
        TMarker* mCx = new TMarker(cxhx[0],cxhx[1],5);
        mCx->SetMarkerColor(4);
        mCx->Draw();
        sprintf(lbl,"Crosshair offset: (%.2f', %.2f')", cxhx[0]-xs[0], cxhx[1]-xs[1]);
        TLatex* L2 = new TLatex(xs[0]-0.6*rr,xs[1]-0.4*rr,lbl);
        L2->SetTextSize(0.03);
        L2->Draw();
        printf("Crosshair at\tH = %s\tV = %s\n",m2dms(cxhx[0]).c_str(),m2dms(cxhx[1]).c_str());
    }
    
    gPad->SetCanvasSize(300,300);
    gPad->Print((fbase+".pdf").c_str());
}

int main(int, char**) {
    string basedir = "/home/mpmendenhall/Documents/aCORN/201501_Alignment";
    vector<string> fnames = listdir(basedir);
    for(auto it = fnames.begin(); it != fnames.end(); it++) {
        if(split(*it," ")[0] != "Finding" || split(*it,".").back() != "csv") continue;
        do_circle_fit(basedir+"/"+dropLast(*it,"."));
    }
}
