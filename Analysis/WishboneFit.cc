#include "WishboneFit.hh"
#include "GraphUtils.hh"
#include "PathUtils.hh"
#include "StringManip.hh"
#include <cassert>
#include <TF2.h>

WishboneFit::WishboneFit(const string& n, OutputManager* pnt): OutputManager(n,pnt) {
    sliceArms[0] = sliceArms[1] = NULL;
}

void WishboneFit::setWishbone(TH2* h) {
    assert(h && !hWishbone);
    hWishbone = h;
    
    // form slices
    for(auto it = hSlices.begin(); it != hSlices.end(); it++) delete *it;
    hSlices = sliceTH2(*hWishbone, X_DIRECTION, true);
    for(int i=0; i<=hWishbone->GetNbinsX()+1; i++)
        hSlices[i]->SetTitle(("Wishbone Slice " + to_str(binE(i)) + " keV").c_str());
    
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
    double E = hWishbone->GetXaxis()->GetBinCenter(bin);
    calcArms(E);
    
    double t0 = 2.8;
    double t1 = 4.5;
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

void WishboneFit::doIt() {
    if(!hWishbone) return;
    
    int b0 = hWishbone->GetXaxis()->FindBin(20);
    int b1 = hWishbone->GetXaxis()->FindBin(750);
    vector<string> hnames;
    for(int b=b0; b<=b1; b++) {
        fitScale(b);
        hSlices[b]->GetXaxis()->SetRangeUser(2,5);
        hSlices[b]->Draw();
        comboFitter.getFitter()->Draw("Same");
        string hname = "Slice_"+to_str(b);
        printCanvas(hname);
        hnames.push_back(plotPath + "/" + hname + ".pdf");
    }
    combo_pdf(hnames,plotPath + "/Slices.pdf");
}

//////////////////
//////////////////

GausWishboneFit::GausWishboneFit(const string& n, OutputManager* pnt):
WishboneFit(n,pnt), MultiGaus(2,"GausWishbone",2), sigma(0.0825) {
    /// initial guess timing
    const size_t npts = 5;
    const double Es[npts] =  {0,   200,  400, 600,  800};
    const double t0s[npts] = {3.3, 3.25, 3.2, 3.2, 3.2};
    const double t1s[npts] = {4.1, 3.95, 3.7, 3.4, 3.2};
    
    for(size_t i=0; i<npts; i++) {
        tau0.SetPoint(i, Es[i], t0s[i]);
        tau1.SetPoint(i, Es[i], t1s[i]);
    }
}

double GausWishboneFit::shapeW(bool upper, double E, double t) const {
    double t0 = upper? tau1.Eval(E) : tau0.Eval(E);
    return exp(-(t-t0)*(t-t0)/(2*sigma*sigma))/(sigma*sqrt(2*M_PI));
}

void GausWishboneFit::fitSliceGaus(int b) {
    double E = binE(b);
    setCenterSigma(0,tau0.Eval(E),sigma);
    setCenterSigma(1,tau1.Eval(E),sigma);
    fit(hSlices[b],false);
}

void GausWishboneFit::doIt() {
    printf("Gaussian wishbone slice fits...\n");
    int nbins = hWishbone->GetXaxis()->FindBin(700);
    TGraphErrors tau0fine(nbins);
    TGraphErrors tau1fine(nbins);
    TGraphErrors sigma0fine(nbins);
    TGraphErrors sigma1fine(nbins);
    for(int b=1; b<=nbins; b++) {
        fitSliceGaus(b);
        double E = binE(b);
        tau0fine.SetPoint(b-1, E, getParameter(1));
        tau0fine.SetPointError(b-1, 0, getParError(1));
        sigma0fine.SetPoint(b-1, E, getParameter(2));
        sigma0fine.SetPointError(b-1, 0, getParError(2));
        tau1fine.SetPoint(b-1, E, getParameter(4));
        tau1fine.SetPointError(b-1, 0, getParError(4));
        sigma1fine.SetPoint(b-1, E, getParameter(5));
        sigma1fine.SetPointError(b-1, 0, getParError(5));
    }
    
    tau1fine.SetTitle("Wishbone timing");
    tau0fine.SetLineColor(4);
    tau1fine.SetLineColor(2);
    sigma0fine.SetLineColor(4);
    sigma1fine.SetLineColor(2);
    
    tau1fine.SetMinimum(3);
    tau1fine.SetMaximum(4.5);
    tau1fine.Draw("AP");
    tau1fine.GetXaxis()->SetTitle("energy [keV]");
    tau1fine.GetYaxis()->SetTitle("proton TOF [#mus]");
    tau1fine.GetYaxis()->SetTitleOffset(1.4);
    tau0fine.Draw("P");
    printCanvas("TauFine");
    
    sigma0fine.SetMinimum(0);
    sigma0fine.SetMaximum(0.2);
    sigma0fine.Draw("AP");
    sigma0fine.SetTitle("Wishbone timing spread");
    sigma0fine.GetXaxis()->SetTitle("energy [keV]");
    sigma0fine.GetYaxis()->SetTitle("proton TOF #sigma [#mus]");
    sigma0fine.GetYaxis()->SetTitleOffset(1.4);
    sigma1fine.Draw("P");
    printCanvas("SigmaFine");
    
    WishboneFit::doIt();
}
