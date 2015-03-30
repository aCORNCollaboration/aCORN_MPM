#ifndef WISHBONEFIT_HH
#define WISHBONEFIT_HH

#include <TH2.h>
#include <TGraph.h>
#include <TMath.h>
#include "LinHistCombo.hh"
#include "OutputManager.hh"
#include "MultiGaus.hh"
#include "SplineFit.hh"
#include <map>
using std::map;

class WeightBins {
public:
    /// constructor
    WeightBins() { }
    /// add content
    void fill(int b, double x, double dx2, double w = 1) { xs[b] += w*x; dx2s[b] += w*w*dx2; ws[b] += w; }
    /// fill content into histogram by bin
    void intoHist(TH1* h);
    map<int,double> xs, dx2s, ws;
};

/// Base class for wishbone plot analysis
class WishboneFit: public OutputManager {
public:
    /// Constructor
    WishboneFit(const string& n, OutputManager* pnt);
    
    /// set up to analyze supplied wishbone plot
    void setWishbone(TH2* h);
    
    /// (unscaled) wishbone upper/lower arm time profile
    virtual double shapeW(bool upper, double E, double t) const = 0;
    /// wishbone model fit in energy bin
    double modelFit(int b, double t) const;
    /// indefinite integral dt at t of unscaled wishbone shape function
    virtual double intShapeW(bool upper, double E, double t) const = 0;
    /// total normalization integral over shape
    virtual double intShapeW(bool upper, double E) const = 0;
    /// characteristic time for specified arm events
    virtual double armTime(bool upper, double E) const = 0;
    /// estimated misidentified tail fraction balance in bin for time cut t
    double tailBalance(int b, double t) const;
    
    /// fit scale factors to wishbone slice; return sum chi^2
    double fitScale(int bin);
    
    /// fit slices with model
    virtual void fitModel();
    /// generate extracted asymmetry plots
    virtual void extractAsymmetry();
    /// inter-arm fill histograms
    virtual void calcGapFill();
    
    int currentBin = 0;                 ///< bin in use for internal analysis
    
protected:
    
    /// calculate arms at given energy
    void calcArms(double E);
    /// energy for wishbone bin
    double binE(int b) const { return hWishbone->GetXaxis()->GetBinCenter(b); }
    /// calculate fit range for energy
    virtual void getComboFitRange(double E, double& t0, double& t1) const;
    /// calculate t0 equal-tails crossover point for energy
    double calcCrossover(int b);
    
    LinHistCombo comboFitter;           ///< fitter for slice scales
    TH2* hWishbone = NULL;              ///< wishbone being analyzed
    vector<TH1F*> hSlices;              ///< wishbone slices for each energy bin
    TH1* sliceArms[2];                  ///< profiles for upper/lower arm slices
    vector<double> C[2];                ///< wishbone arm magnitudes for each slice
    vector<double> dC[2];               ///< fit uncertainty on arm magnitudes
    TGraph gt0;                         ///< optimal t0 cut point
};

/// Wishbone fit assuming Gaussian profile
class GausWishboneFit: public WishboneFit, protected MultiGaus {
public:
    /// Constructor
    GausWishboneFit(const string& n, OutputManager* pnt);
    
    /// (unscaled) wishbone upper/lower arm time profile
    virtual double shapeW(bool upper, double E, double t) const;
    /// indefinite integral dt at t of unscaled wishbone shape function
    virtual double intShapeW(bool uppter, double E, double t) const;
    /// total normalization integral over shape
    virtual double intShapeW(bool, double) const { return 1.; }
    /// characteristic time for specified arm events
    virtual double armTime(bool upper, double E) const { return tau[upper].mySpline.Eval(E); }
    
    /// evaluation for fit of single slice
    double sliceFitEval(double* x, double* p);
    
    /// perform overall fitting / asymmetry extraction routine
    virtual void fitModel();
    
    SplineFit tau[2];   ///< wishbone arm central timing [us]
    SplineFit sigma[2]; ///< wishbone arm timing width [us]
protected:
    /// 2-peak gaussian fit to a single slice
    void fitSliceGaus(int b);
    /// calculate fit range for energy
    virtual void getComboFitRange(double E, double& t0, double& t1) const;
};

#endif
