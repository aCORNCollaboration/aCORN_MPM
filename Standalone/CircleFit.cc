#include "StringManip.hh"
#include "GraphicsUtils.hh"
#include "Matrix.hh"

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
            double c = xs[i]-params[0];
            double s = ys[i]-params[1];
            double r0 = sqrt(c*c+s*s);
            c /= r0;
            s /= r0;
            double r = 1./sqrt((iSigma(0,0)*c+iSigma(0,1)*s)*c + (iSigma(1,0)*c+iSigma(1,1)*s)*s);
            s_err += pow(r0-r,2);
        }
        return s_err;
    }
    
    vector<double> xs;
    vector<double> ys;
    Matrix<2,2,double> iSigma; ///< inverse covariance matrix
    
    /// calculate initial guess
    void initGuess(double& x0, double& y0, double& r0) {
        x0 = y0 = 0;
        for(size_t i=0; i<xs.size(); i++) {
            x0 += xs[i];
            y0 += ys[i];
        }
        x0 /= xs.size();
        y0 /= ys.size();
        double params[5] = {x0, y0, 1e-6, 0, 1e-6};
        r0 = sqrt(circleMin(params)/xs.size());
    }
    
    TGraph* ptsGraph() const {
        TGraph* g = new TGraph(xs.size());
        for(size_t i=0; i<xs.size(); i++)
            g->SetPoint(i,xs[i],ys[i]);
        return g;
    }
};

int main(int, char**) {
    CircleMin QM;
    std::ifstream t("/home/mpmendenhall/Downloads/Finding the Center of the Thick Teflon Bottom of EM_12.30.2014.csv");
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
    double rms = QM.circleMin(xs);
    cout << "Minimum: f( ";
    for(unsigned int i=0; i<nvar; i++) cout << xs[i] << " ";
    cout << "): rms = " << rms << "\n";
    
    
    TGraph gCirc(100);
    double rr = sqrt(0.5*(xs[2]+xs[4]));
    for(int i=0; i<100; i++) {
        double th = i*2*M_PI/99.;
        gCirc.SetPoint(i, xs[0]+rr*cos(th), xs[1]+rr*sin(th));
    }
    gCirc.SetLineColor(4);
    gCirc.Draw("AC");
    gCirc.SetTitle("theodolite circle fit");
    gCirc.GetXaxis()->SetTitle("horizontal [arcmin]");
    gCirc.GetYaxis()->SetTitle("vertical [arcmin]");
    gCirc.GetYaxis()->SetTitleOffset(1.5);
        
    TPolyLine* LE = makeEllipse(xs[0], xs[1], &QM.iSigma[0]);
    LE->SetLineColor(2);
    LE->Draw();
    
    TGraph* g = QM.ptsGraph();
    g->Draw("*");
    
    TMarker mCenter(xs[0],xs[1],2);
    TMarker mCx(cxhx[0],cxhx[1],5);
    mCenter.SetMarkerColor(2);
    mCenter.Draw();
    mCx.SetMarkerColor(4);
    mCx.Draw();
    
    char lbl[1024];
    sprintf(lbl,"#splitline{   Center (%.2f,%.2f), RMS %.2f}{r_{xx} = %.2f; r_{yy} = %.2f; r_{xy} = %.2f}", xs[0], xs[1], rms, sqrt(xs[2]), sqrt(xs[4]), sqrt(xs[3]));
    TLatex L(xs[0]-0.7*rr,xs[1]+40,lbl);
    L.SetTextSize(0.03);
    L.Draw();
    
    sprintf(lbl,"Crosshair: (%.2f,%.2f)", cxhx[0], cxhx[1]);
    TLatex L2(xs[0]-0.6*rr,xs[1]-40,lbl);
    L2.SetTextSize(0.03);
    L2.Draw();
    
    gPad->SetCanvasSize(300,300);
    gPad->Print("circle.pdf");
    return 0;
}