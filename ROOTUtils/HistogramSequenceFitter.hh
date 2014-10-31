#ifndef HISTOGRAMSEQUENCEFITTER_HH
#define HISTOGRAMSEQUENCEFITTER_HH

#include <cmath>
#include <vector>

#include <TH1.h>
#include <TF1.h>
#include <TF2.h>
#include <TFitResultPtr.h>
#include <TGraph2DErrors.h>
#include <Math/GSLMinimizer.h>
#include <Math/Functor.h>

using std::vector;
using std::pair;
typedef vector< pair<double,double> > intervalList;

/// base class for fitting function given integrals over intervals
class IntervalIntegralFitter {
public:
    /// Constructor
    IntervalIntegralFitter(unsigned int npar);
    /// Destructor
    virtual ~IntervalIntegralFitter();
    
    /// calculate fit results
    void fit();
    
    const unsigned int N;               ///< number of parameters
    double* myParams;                   ///< parameter values
    double* mySteps;                    ///< optimization step sizes for parameters
    
    vector<intervalList> intervals;     ///< intervals over which function is evaluated
    vector<double> integrals;           ///< integral value over each interval
    vector<double> dIntegrals;          ///< uncertainties on integral values
    
    /// sum-squared-error between fit and values
    double eval_error(const double* params) const;
    
    /// definite integral over intervals (using supplied or internal parameter list)
    double operator()(const intervalList& L, const double* params = NULL) const;
    
    /// indefinite integral of fit function (subclass me!)
    virtual double integ_f(double t, const double* params) const = 0;
    
    /// initial fit parameter,step setup
    virtual void init_params() = 0;
    
protected:
    
    ROOT::Math::GSLMinimizer myMin;
    ROOT::Math::Functor myf;
};

/// Exponential interval fitter
class ExponentialIntegralFitter: public IntervalIntegralFitter {
public:
    /// constructor
    ExponentialIntegralFitter(): IntervalIntegralFitter(2), T0(0) { }
    
    /// indefinite integral of variation function (subclass me!)
    virtual double integ_f(double t, const double* params) const;
    
    /// initial fit parameter estimation
    virtual void init_params();
    
protected:
    
    static TF1 expFit;  ///< fitter for parameter estimation
    double T0;          ///< fit reference start time
};

/// Polynomially varying sequence
class PolynomialIntegralFitter: public IntervalIntegralFitter {
public:
    /// constructor
    PolynomialIntegralFitter(unsigned int npar): IntervalIntegralFitter(npar), T0(0) { }
    /// indefinite integral of variation function (subclass me!)
    virtual double integ_f(double t, const double* params) const;
    
    /// initial fit parameter estimation
    virtual void init_params();
    
protected:
    
    double T0;          ///< fit reference start time
};


/////////////////////////////////////////
/////////////////////////////////////////
/////////////////////////////////////////

/// Bin-by-bin fits for histograms integrating over a time-evolving quantity
class HistogramSequenceFitter {
public:
    /// constructor
    HistogramSequenceFitter(IntervalIntegralFitter* f): myFitter(f) { }
    /// destructor
    virtual ~HistogramSequenceFitter() {}
    
    /// add fit data point histogram
    void addData(const TH1* h, const intervalList& dt);
    /// calculate fit results
    void fit();
    /// generate interpolated histogram for time range
    TH1* interpolate(const intervalList& dt, TH1* h = NULL) const;
    
protected:
    
    IntervalIntegralFitter* myFitter;   ///< fitter function operating on this class
    vector<const TH1*> hs;              ///< input data point histograms
    vector<intervalList> dts;           ///< time intervals
    
    vector< vector<double> > fts;       ///< fit parameters for each bin
};


#endif
