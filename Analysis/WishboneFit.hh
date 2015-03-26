#ifndef WISHBONEFIT_HH
#define WISHBONEFIT_HH

#include <TH2.h>
#include <TGraph.h>
#include <TMath.h>
#include "LinHistCombo.hh"
#include "OutputManager.hh"
#include "MultiGaus.hh"
#include "SplineFit.hh"

/// Base class for wishbone plot analysis
class WishboneFit: public OutputManager {
public:
    /// Constructor
    WishboneFit(const string& n, OutputManager* pnt);
    
    /// set up to analyze supplied wishbone plot
    void setWishbone(TH2* h);
    
    /// (unscaled) wishbone upper/lower arm time profile
    virtual double shapeW(bool upper, double E, double t) const = 0;
    
    /// fit scale factors to wishbone slice; return sum chi^2
    double fitScale(int bin);
    
    /// perform overall fitting / asymmetry extraction routine
    virtual void doIt();
    
protected:
    
    /// calculate arms at given energy
    void calcArms(double E);
    /// energy for wishbone bin
    double binE(int b) const { return hWishbone->GetXaxis()->GetBinCenter(b); }
    
    LinHistCombo comboFitter;           ///< fitter for slice scales
    TH2* hWishbone = NULL;              ///< wishbone being analyzed
    vector<TH1F*> hSlices;              ///< wishbone slices for each energy bin
    TH1* sliceArms[2];                  ///< profiles for upper/lower arm slices
    vector<double> C[2];                ///< wishbone arm magnitudes for each slice
};

/// Wishbone fit assuming Gaussian profile
class GausWishboneFit: public WishboneFit, protected MultiGaus {
public:
    /// Constructor
    GausWishboneFit(const string& n, OutputManager* pnt);
    
    /// (unscaled) wishbone upper/lower arm time profile
    virtual double shapeW(bool upper, double E, double t) const;
    
    /// evaluation for fit of single slice
    double sliceFitEval(double* x, double* p);
    
    /// perform overall fitting / asymmetry extraction routine
    virtual void doIt();
    
    SplineFit tau[2];   ///< wishbone arm central timing [us]
    SplineFit sigma[2]; ///< wishbone arm timing width [us]
protected:
    /// 2-peak gaussian fit to a single slice
    void fitSliceGaus(int b);
};

#endif
