#ifndef MULTIGAUS_HH
#define MULTIGAUS_HH 1

#include <string>
using std::string;
#include <vector>
using std::vector;

#include <climits>

#include <TF1.h>
#include <TH1F.h>
#include "FloatErr.hh"

/// class for fitting multi-peak gaussians
class MultiGaus {
public:
    
    /// correlated subpeaks specification
    struct corrPeak {
        unsigned int mainPeak;
        double relCenter;
        double relHeight;
        double relWidth;
    };
    
    /// constructor
    MultiGaus(unsigned int n, const string& name, float ns = 1.5): nSigma(ns), npks(n), iguess(new double[3*n]), myTF1(new TF1(name.c_str(),this,0,0,3*n)) { }
    
    /// destructor
    ~MultiGaus();
    
    /// add correlated peak
    void addCorrelated(unsigned int n, double relCenter, double relHeight, double relWidth = 0);
    /// fill initial values array
    void setParameter(unsigned int n, double p);
    /// set center,sigma initial values
    void setCenterSigma(unsigned int n, double c, double s);
        
    /// get fit parameter
    double getParameter(unsigned int n) const;
    /// get fit parameter error
    double getParError(unsigned int n) const;
    /// get parameter+error as float_err
    float_err getPar(unsigned int n) const;
    /// display fit results
    void display() const;
    
    /// get TF1 with appropriate pre-set values
    TF1* getFitter();
    
    /// single-peak fit estimate, assuming center and sigma initial guess set
    void fitEstimate(TH1* h, unsigned int n = INT_MAX);
    /// fit a TH1 after initial peak centers/widths have been guessed; update inital guess
    void fit(TH1* h, bool draw = true);
    
    /// gaussian evaluation function
    double operator() (double* x, double* par);
    
    float nSigma;               ///< number of sigma peak width to fit
    
    const unsigned int npks;    ///< number of peaks being fitted
    
protected:
    double* iguess;             ///< inital guess at peak positions
    TF1* myTF1;                 ///< TF1 using this class as its fit function
    vector<corrPeak> corrPeaks; ///< correlated subpeaks
};

int iterGaus(TH1* h0, TF1* gf, unsigned int nit, float mu, float sigma, float nsigma = 1.5, float asym = 0);

#endif
