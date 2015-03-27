#include "WishboneFit.hh"
#include "GraphUtils.hh"
#include "GraphicsUtils.hh"
#include "PathUtils.hh"
#include "StringManip.hh"
#include <cassert>
#include <TF2.h>
#include <TMath.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

WishboneFit::WishboneFit(const string& n, OutputManager* pnt): OutputManager(n,pnt) {
    sliceArms[0] = sliceArms[1] = NULL;
}

void WishboneFit::setWishbone(TH2* h) {
    assert(h && !hWishbone);
    hWishbone = h;
    
    // form slices
    for(auto it = hSlices.begin(); it != hSlices.end(); it++) delete *it;
    hSlices = sliceTH2(*hWishbone, X_DIRECTION, true);
    for(int i=0; i<=hWishbone->GetNbinsX()+1; i++) {
        hSlices[i]->SetTitle(("Wishbone Slice " + to_str(binE(i)) + " keV").c_str());
        hSlices[i]->GetYaxis()->SetTitle("rate [Hz/MeV/#mus]");
        hSlices[i]->GetYaxis()->SetTitleOffset(1.35);
    }
    
    // set up arm profiles
    // TODO comboFitter.Clear();
    for(int i=0; i<2; i++) {
        if(sliceArms[i]) delete sliceArms[i];
        sliceArms[i] = (TH1*)hSlices[0]->Clone(("sliceArms_"+to_str(i)).c_str());
        comboFitter.addTerm(sliceArms[i]);
        C[i].clear();
        C[i].resize(hSlices.size());
    }
}

void WishboneFit::calcArms(double E) {
    for(int i=0; i<2; i++)
        for(int b=1; b<=sliceArms[i]->GetNbinsX(); b++)
            sliceArms[i]->SetBinContent(b, shapeW(i,E,sliceArms[i]->GetBinCenter(b)));
}

double WishboneFit::fitScale(int bin) {
    assert(hWishbone);
    double E = binE(bin);
    calcArms(E);
    
    double t0, t1;
    getComboFitRange(E,t0,t1);
    comboFitter.Fit(hSlices[bin], t0, t1);
    for(int i=0; i<2; i++) C[i][bin] = comboFitter.coeffs[i];
    
    // calculate residuals chi^2
    const TAxis* AX = hSlices[bin]->GetXaxis();
    int b0 = AX->FindBin(t0);
    int b1 = AX->FindBin(t1);
    double chisq = 0;
    for(int b=b0; b<=b1; b++) {
        double dy = (hSlices[bin]->GetBinContent(b) - comboFitter.eval(AX->GetBinCenter(b)))/hSlices[bin]->GetBinError(b);
        chisq += dy*dy;
    }
    return chisq;
}

void WishboneFit::fitModel() {
    if(!hWishbone) return;
    
    int b0 = hWishbone->GetXaxis()->FindBin(20);
    int b1 = hWishbone->GetXaxis()->FindBin(750);
    vector<string> hnames;
    for(int b=b0; b<=b1; b++) {
        fitScale(b);
        double cx = calcCrossover(b);
        gt0.SetPoint(b,binE(b),cx);
                     
        hSlices[b]->GetXaxis()->SetRangeUser(2,5);
        hSlices[b]->SetMinimum(-0.2);
        if(hSlices[b]->GetMaximum() < 1.2) hSlices[b]->SetMaximum(1.2);
        hSlices[b]->Draw();
        drawVLine(cx,defaultCanvas,6);
        comboFitter.getFitter()->Draw("Same");
        string hname = "Slice_"+to_str(b);
        printCanvas(hname);
        hnames.push_back(plotPath + "/" + hname + ".pdf");
    }
    combo_pdf(hnames,plotPath + "/Slices.pdf");
}

void WishboneFit::getComboFitRange(double, double& t0, double& t1) const {
    t0 = 2.8;
    t1 = 4.5;
}

double WishboneFit::tailBalance(int b, double t) const {
    double E = binE(b);
    double tails = -C[0][b]*(intShapeW(0, E, 10) - intShapeW(0, E, t));
    tails += C[1][b]*(intShapeW(1, E, t) - intShapeW(1, E, 0));
    return tails;
}

double wishbone_tails(double t, void* params) {
    WishboneFit* WBF = (WishboneFit*)params;
    return WBF->tailBalance(WBF->currentBin, t);
}

double WishboneFit::calcCrossover(int b) {
    currentBin = b;
    double E = binE(b);
    if(C[0][b] <= 0 || C[1][b] <= 0) return 0;
    
    double r = 0;
    double x_lo, x_hi;
    getComboFitRange(E, x_lo, x_hi);
    
    gsl_function F;
    F.function = &wishbone_tails;
    F.params = this;
    gsl_root_fsolver* s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
        
    int iter = 0, max_iter = 100;
    int status;
    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        r = gsl_root_fsolver_root(s);
        x_lo = gsl_root_fsolver_x_lower(s);
        x_hi = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(x_lo, x_hi, 0, 0.0001);
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(s);
    
    return r;
}

//////////////////
//////////////////

GausWishboneFit::GausWishboneFit(const string& n, OutputManager* pnt):
WishboneFit(n,pnt), MultiGaus(2,"GausWishbone",2) {
    /// initial guess timing
    const size_t npts = 12;
    const double Es[npts] =  {20,   40,   70,   100,  200,  300,  400,  500,  550, 600, 650,  700};
    const double t0s[npts] = {3.3,  3.3,  3.3,  3.25, 3.25, 3.2,  3.2,  3.2,  3.2, 3.2, 3.2,  3.2};
    const double t1s[npts] = {4.1,  4.0,  4.0,  4.0,  3.95, 3.8,  3.7,  3.5,  3.5, 3.4, 3.4,  3.2};
    
    for(size_t i=0; i<npts; i++) {
        for(int j=0; j<2; j++) {
            tau[j].mySpline.SetPoint(i, Es[i], (j?t1s:t0s)[i]);
            sigma[j].mySpline.SetPoint(i, Es[i], 0.083);
        }
    }
}

double GausWishboneFit::shapeW(bool upper, double E, double t) const {
    double t0 = tau[upper].mySpline.Eval(E);
    double s0 = sigma[upper].mySpline.Eval(E);
    return exp(-(t-t0)*(t-t0)/(2*s0*s0))/(s0*sqrt(2*M_PI));
}

double GausWishboneFit::intShapeW(bool upper, double E, double t) const {
    double t0 = tau[upper].mySpline.Eval(E);
    double s0 = sigma[upper].mySpline.Eval(E);
    return TMath::Erf( (t-t0)/s0 );
}


void GausWishboneFit::fitSliceGaus(int b) {
    double E = binE(b);
    for(int i=0; i<2; i++)
        setCenterSigma(i, tau[i].mySpline.Eval(E), sigma[i].mySpline.Eval(E));
    fit(hSlices[b],false);
}

void GausWishboneFit::fitModel() {
    printf("Gaussian wishbone slice fits...\n");
    int nbins = hWishbone->GetXaxis()->FindBin(700);
    TGraphErrors taufine[2];
    TGraphErrors sigmafine[2];
    for(int b=1; b<=nbins; b++) {
        fitSliceGaus(b);
        double E = binE(b);
        for(int i=0; i<2; i++) {
            taufine[i].SetPoint(b-1, E, getParameter(1+3*i));
            taufine[i].SetPointError(b-1, 0, getParError(1+3*i));
            sigmafine[i].SetPoint(b-1, E, getParameter(2+3*i));
            sigmafine[i].SetPointError(b-1, 0, getParError(2+3*i));
        }
    }
    
    // spline fit timing adjustment
    for(int i=0; i<2; i++) {
        TF1* fTau = tau[i].getFitter(&taufine[i]);
        taufine[i].Fit(fTau,"RN");
        tau[i].updateSpline();
        
        TF1* fSigma = sigma[i].getFitter(&sigmafine[i]);
        sigmafine[i].Fit(fSigma,"RN");
        sigma[i].updateSpline();
    }
    
    WishboneFit::fitModel();
    
    // plot setup
    for(int i=0; i<2; i++) {
        taufine[i].SetTitle("Wishbone timing");
        taufine[i].SetLineColor(4-2*i);
        sigmafine[i].SetLineColor(4-2*i);
        taufine[i].SetMinimum(3);
        taufine[i].SetMaximum(4.5);
        sigmafine[i].SetMinimum(0);
        sigmafine[i].SetMaximum(0.2);
    }
    taufine[1].Draw("AP");
    taufine[1].GetXaxis()->SetTitle("energy [keV]");
    taufine[1].GetYaxis()->SetTitle("proton TOF [#mus]");
    taufine[1].GetYaxis()->SetTitleOffset(1.4);
    taufine[0].Draw("P");
    gt0.SetMarkerStyle(6);
    gt0.Draw("P");
    //
    for(int i=0; i<2; i++) {
        tau[i].mySpline.SetLineColor(6);
        tau[i].mySpline.SetLineWidth(2);
        tau[i].mySpline.Draw("PC");
    }
    //
    printCanvas("TauFine");
    
    sigmafine[0].Draw("AP");
    sigmafine[0].SetTitle("Wishbone timing spread");
    sigmafine[0].GetXaxis()->SetTitle("energy [keV]");
    sigmafine[0].GetYaxis()->SetTitle("proton TOF #sigma [#mus]");
    sigmafine[0].GetYaxis()->SetTitleOffset(1.4);
    sigmafine[1].Draw("P");
    //
    for(int i=0; i<2; i++) {
        sigma[i].mySpline.SetLineColor(6);
        sigma[i].mySpline.SetLineWidth(2);
        sigma[i].mySpline.Draw("PC");
    }
    //
    printCanvas("SigmaFine");
    
    
}

void GausWishboneFit::getComboFitRange(double E, double& t0, double& t1) const {
    t0 = 2.8;
    t1 = tau[1].mySpline.Eval(E)+0.5;
}
