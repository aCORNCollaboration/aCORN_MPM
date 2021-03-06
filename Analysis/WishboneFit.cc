/// \file WishboneFit.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

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

void WeightBins::intoHist(TH1* h) {
    for(auto const& bv: xs) h->SetBinContent(bv.first, bv.second);
    for(auto const& bv: dx2s) h->SetBinError(bv.first, sqrt(bv.second));
}

/////////////////////
/////////////////////
/////////////////////

void WishboneSeparator::setWishbone(TH2* h) {
    assert(h && !hWishbone);
    hWishbone = h;
    
    for(auto sl: hSlices) delete sl;
    hSlices = sliceTH2(*hWishbone, X_DIRECTION, true);
    for(int i=0; i<=hWishbone->GetNbinsX()+1; i++) {
        hSlices[i]->SetTitle(("Wishbone Slice " + to_str(binE(i)) + " keV").c_str());
    }
}

WishboneSeparator::~WishboneSeparator() { 
    for(auto sl: hSlices) delete sl;
}

void WishboneSeparator::extractAsymmetry() {
    int npts = 0;
    int nptsa = 0;
    for(int b=1; b<(int)hSlices.size(); b++) {
        double E = binE(b);
        double tm = gt0.Eval(E);

        double t0,t1;
        getComboFitRange(E, t0, t1);
        double dadm, dadp;
        double adm = integralAndErrorInterp(hSlices[b],t0,tm,dadm,true);
        double adp = integralAndErrorInterp(hSlices[b],tm,t1,dadp,true);
        dataN[0].SetPoint(npts, E, adm);
        dataN[1].SetPoint(npts, E, adp);
        dataN[0].SetPointError(npts, 0, dadm);
        dataN[1].SetPointError(npts, 0, dadp);
        if(adm+adp) {
            double ad = (adp-adm)/(adp+adm);
            double dad = sqrt(dadm*dadm + dadp*dadp)/(adm+adp);
            dataA.SetPoint(nptsa, E, 100*ad);
            dataA.SetPointError(nptsa, 0, 100*dad);
            nptsa++;
        }
        npts++;
    }
    
    for(int i=0; i<2; i++) {
        dataN[i].SetMarkerStyle(24+3*i);
        dataN[i].SetMarkerSize(0.2);
        dataN[i].SetMarkerColor(2);
        dataN[i].SetLineColor(2);
    }
    
    dataA.SetMarkerStyle(7);
    dataA.SetMarkerColor(2);
    dataA.SetLineColor(2);
}

void WishboneSeparator::compareAsym(WishboneSeparator& WB) {
    WB.dataA.SetMarkerColor(4);
    WB.dataA.SetLineColor(4);
    
    TF1 lineFit("lineFit","pol0",100,400);
    WB.dataA.Fit(&lineFit,"RN");
    dataA.Fit(&lineFit,"RN");
    
    dataA.Draw("AP");
    WB.dataA.Draw("P");
    printCanvas("AsymCompare");
    
    gt0.SetLineColor(2);
    WB.gt0.SetLineColor(4);
    WB.gt0.SetMarkerColor(4);
    WB.gt0.SetMarkerStyle(7);
    gt0.Draw("AC");
    gt0.SetTitle("Wishbone timing cuts comparison");
    gt0.SetMinimum(3.2);
    gt0.SetMaximum(3.8);
    gt0.GetXaxis()->SetTitle("energy [keV]");
    gt0.GetXaxis()->SetRangeUser(0,600);
    gt0.GetYaxis()->SetTitle("timing cut [#mus]");
    WB.gt0.Draw("CP");
    printCanvas("T0Compare");
}

double WishboneSeparator::fitAsym(double E0, double E1, double& err) {
    TF1 fLin("fLin","pol0",E0,E1);
    dataA.Fit(&fLin, "RN");
    err = fLin.GetParError(0);
    return fLin.GetParameter(0);
}

/////////////////////
/////////////////////
/////////////////////

ManualWishboneSeparator::ManualWishboneSeparator(const string& n, OutputManager* pnt): WishboneSeparator(n,pnt) {
    const size_t npts = 7;
    const double E0s[npts] = {0,    100,  200,  300,  400,  500,  600 };
    const double t0s[npts] = {3.72, 3.64, 3.55, 3.48, 3.40, 3.31, 3.25 };
    for(size_t i=0; i<npts; i++) gt0.SetPoint(i, E0s[i], t0s[i]);
}
 
void ManualWishboneSeparator::getComboFitRange(double E, double& t0, double& t1) const {
    t0 = 3.05;
    t1 = 4.3 - 0.8*(E/700)*(E/700);
}
 
/////////////////////
/////////////////////
/////////////////////

WishboneMidcentroid::WishboneMidcentroid(const string& n, OutputManager* pnt): WishboneSeparator(n,pnt) {
    const size_t npts = 21;
    const double t0s[npts] = {3.65, 3.6375, 3.6, 3.6125, 3.6, 3.5875, 3.5625, 3.55, 3.55, 3.5375, 3.5125, 3.5, 3.5, 3.4875, 3.475, 3.45, 3.4375, 3.425, 3.4125, 3.3875, 3.375};
    for(size_t i=0; i<npts; i++) gt0.SetPoint(i, 110+20*i, t0s[i]);
}

/////////////////////
/////////////////////
/////////////////////

WishboneFit::WishboneFit(const string& n, OutputManager* pnt): WishboneSeparator(n,pnt) {
    sliceArms[0] = sliceArms[1] = nullptr;
}

void WishboneFit::setWishbone(TH2* h) {
    WishboneSeparator::setWishbone(h);

    // set up arm profiles
    // TODO comboFitter.Clear();
    for(int i=0; i<2; i++) {
        if(sliceArms[i]) delete sliceArms[i];
        sliceArms[i] = (TH1*)hSlices[0]->Clone(("sliceArms_"+to_str(i)).c_str());
        comboFitter.addTerm(sliceArms[i]);
        C[i].clear();
        C[i].resize(hSlices.size());
        dC[i].clear();
        dC[i].resize(hSlices.size());
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
    for(int i=0; i<2; i++) {
        C[i][bin] = comboFitter.coeffs[i];
        dC[i][bin] = comboFitter.dcoeffs[i];
    }
    
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

double WishboneFit::modelFit(int b, double t) const {
    if(b < 0 || b >= (int)hSlices.size()) return 0;
    double E = binE(b);
    return C[0][b]*shapeW(0,E,t) + C[1][b]*shapeW(1,E,t);
}

void WishboneFit::fitModel() {
    if(!hWishbone) return;
    
    int b0 = hWishbone->GetXaxis()->FindBin(20);
    int b1 = hWishbone->GetXaxis()->FindBin(750);
    TF1* cf = comboFitter.getFitter();
    cf->SetNpx(300);
    defaultCanvas.SetRightMargin(0.04);
    
    vector<string> hnames;
    for(int b=b0; b<=b1; b++) {
        fitScale(b);
        double cx = calcCrossover(b);
        gt0.SetPoint(b,binE(b),cx);
                     
        hSlices[b]->GetXaxis()->SetRangeUser(2,5);
        hSlices[b]->SetMinimum(-0.2);
        if(hSlices[b]->GetMaximum() < 1.2) hSlices[b]->SetMaximum(1.2);
        hSlices[b]->Draw();
        addDeletable(drawVLine(cx, &defaultCanvas,6));
        cf->Draw("Same");
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
    if(tailBalance(b,x_lo) >= 0 || tailBalance(b,x_hi) <= 0) return 0;
       
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

void WishboneFit::extractAsymmetry() {
    
    WishboneSeparator::extractAsymmetry();
    
    TGraphErrors modelN[2];
    TGraphErrors modelA;
    TGraph dtSens;
    int npts = 0;
    TF1 lineFit("lineFit","pol0",100,400);
    
    for(int b=1; b<(int)hSlices.size(); b++) {
        if(C[0][b] <= 0 || C[1][b] <= 0) continue;
        double E = binE(b);
        
        double tm = gt0.Eval(E);
        if(!(armTime(0,E) < tm && tm < armTime(1,E))) continue;
        
        double n0 = intShapeW(0, E);
        double n1 = intShapeW(1, E);
        double amm = C[0][b]*n0;
        double amp = C[1][b]*n1;
        modelN[0].SetPoint(npts, E, amm);
        modelN[1].SetPoint(npts, E, amp);
        modelN[0].SetPointError(npts, 0, dC[0][b]*n0);
        modelN[1].SetPointError(npts, 0, dC[1][b]*n0);
        double am = (amp-amm)/(amp+amm);
        double dam = sqrt(dC[0][b]*dC[0][b]*n0*n0 + dC[1][b]*dC[1][b]*n1*n1)/(amm+amp);
        modelA.SetPoint(npts, E, 100*am);
        modelA.SetPointError(npts, 0, 100*dam);
        
        double sens = 2*modelFit(b,tm)/(amm+amp);
        dtSens.SetPoint(npts, E, 100*sens/20);
        
        npts++;
    }
    
    for(int i=0; i<2; i++) {
        modelN[i].SetMarkerStyle(24+3*i);
        modelN[i].SetMarkerSize(0.2);
        modelN[i].SetMarkerColor(4);
        modelN[i].SetLineColor(4);
    }
    modelA.SetMarkerColor(4);
    modelA.SetLineColor(4);
    
    defaultCanvas.SetRightMargin(0.04);
    
    modelN[0].Draw("AP");
    modelN[0].SetTitle("aCORN wishbone counts");
    modelN[0].GetXaxis()->SetTitle("energy [keV]");
    modelN[0].GetXaxis()->SetRangeUser(0,600);
    modelN[0].GetYaxis()->SetTitle("rate [Hz/MeV]");
    modelN[1].Draw("P");
    dataN[0].Draw("P");
    dataN[1].Draw("P");
    printCanvas("WishboneCounts");
    
    dataA.SetMinimum(0);
    dataA.SetMaximum(20);
    dataA.SetMarkerColor(2);
    dataA.SetLineColor(2);
    dataA.SetMarkerStyle(7);
    dataA.Draw("AP");
    dataA.SetTitle("aCORN wishbone asymmetry");
    dataA.GetXaxis()->SetTitle("energy [keV]");
    dataA.GetXaxis()->SetRangeUser(0,500);
    dataA.GetYaxis()->SetTitle("branch asymmetry [%]");
    modelA.Draw("P");
    
    dtSens.SetMarkerStyle(7);
    //dtSens.Draw("PL");
    
    printCanvas("Asymmetry");
    
    modelA.Fit(&lineFit,"RN");
    dataA.Fit(&lineFit,"RN");
}

void WishboneFit::calcGapFill() {
    WeightBins wbx, wby;
    for(int i=0; i<hWishbone->GetNcells(); i++) {
        Int_t binx, biny, binz;
        hWishbone->GetBinXYZ(i,binx,biny,binz);
        if(binx == 0 || binx > hWishbone->GetNbinsX() || biny == 0 || biny > hWishbone->GetNbinsY()) continue;
        double E = binE(binx);
        double dE = hWishbone->GetXaxis()->GetBinWidth(binx);
        double t = hWishbone->GetYaxis()->GetBinCenter(biny);
        double dt = hWishbone->GetYaxis()->GetBinWidth(biny);
        if(t < armTime(false, E) || t > armTime(true, E)) continue;
        double zexp = modelFit(binx,t);
        if(zexp > 0.05) continue;
        double z = hWishbone->GetBinContent(i) - zexp;
        double dz = hWishbone->GetBinError(i);
        wbx.fill(binx, z, dz*dz, 1000*dt);
        if(E>100) wby.fill(biny, z, dz*dz, dE);
    }
    
    TH1F* fillx = axisHist(*hWishbone, "hWishboneFill_E", X_DIRECTION);
    fillx->SetTitle("Wishbone inter-arm fill");
    fillx->GetYaxis()->SetTitle("rate [mHz/MeV]");
    fillx->GetXaxis()->SetRangeUser(0,600);
    fillx->SetMaximum(50);
    
    TH1F* filly = axisHist(*hWishbone, "hWishboneFill_t", Y_DIRECTION);
    filly->GetYaxis()->SetTitle("rate [mHz/#mus]");
    filly->SetTitle("Wishbone inter-arm fill");
    filly->GetXaxis()->SetRangeUser(3,4);
    wbx.intoHist(fillx);
    wby.intoHist(filly);
    
    fillx->Draw();
    addDeletable(drawHLine(0, &defaultCanvas, 1, 2));
    printCanvas("hWishboneFill_E");
    
    filly->Draw();
    addDeletable(drawHLine(0, &defaultCanvas, 1, 2));
    printCanvas("hWishboneFill_t");
    
    delete fillx;
    delete filly;
}

//////////////////
//////////////////

GausWishboneFit::GausWishboneFit(const string& n, OutputManager* pnt):
WishboneFit(n,pnt), MultiGaus(2,"GausWishbone",2) {
    /// initial guess timing
    const size_t npts = 11;
    const double Es[npts] =  {20,   40,   70,   100,  200,  300,  400,  500,  550, 600, 640};
    const double t0s[npts] = {3.3,  3.3,  3.3,  3.25, 3.25, 3.2,  3.2,  3.2,  3.2, 3.2, 3.2};
    const double t1s[npts] = {4.1,  4.0,  4.0,  4.0,  3.95, 3.8,  3.7,  3.5,  3.5, 3.4, 3.3};
    
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
    return (1+TMath::Erf( (t-t0)/s0 ))*0.5;
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
    defaultCanvas.SetRightMargin(0.04);
    for(int i=0; i<2; i++) {
        taufine[i].SetTitle("Wishbone timing");
        taufine[i].SetLineColor(4-2*i);
        sigmafine[i].SetLineColor(4-2*i);
        taufine[i].SetMinimum(3);
        taufine[i].SetMaximum(4.5);
        sigmafine[i].SetMinimum(0);
        sigmafine[i].SetMaximum(0.2);
    }
    defaultCanvas.SetRightMargin(0.15);
    makeRBpalette();
    hWishbone->Draw("Col Z");
    taufine[1].Draw("P");
    taufine[1].GetXaxis()->SetTitle("energy [keV]");
    taufine[1].GetYaxis()->SetTitle("proton TOF [#mus]");
    taufine[0].Draw("P");
    gt0.SetMarkerStyle(8);
    gt0.SetMarkerSize(0.2);
    gt0.SetMarkerColor(4);
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
