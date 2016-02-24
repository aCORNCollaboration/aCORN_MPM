/// \file CircleMin.hh Fit circle (ellipse) to set of points
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2016
 
#include "Matrix.hh"
#include <TGraph.h>
#include <Math/GSLMinimizer.h>

/// Ellipse fit minimization function
class CircleMin {
public:
    CircleMin(): min( ROOT::Math::kVectorBFGS ) { }
    
    /// minimization function, params: x,y,rxx,rxy,ryy
    double circleMin(const double* params);
    
    /// add point to fit data
    void addPoint(double x, double y) { xs.push_back(x); ys.push_back(y); }
    
    bool verbose = false;               ///< verbose printout on evaluation step
    Matrix<2,2,double> iSigma;          ///< inverse covariance matrix
    ROOT::Math::GSLMinimizer min;       ///< the minimizer
    
    /// transform points p -> M*(p-p0)
    void transform(double x0, double y0, const Matrix<2,2,double>& M);
    
    /// perform fit
    double doFit();
    /// calculate initial guess
    void initGuess(double& x0, double& y0, double& r0);
    /// produce TGraph showing points
    TGraph* ptsGraph() const;
    
    vector<double> xs;          ///< x coordinate each point
    vector<double> ys;          ///< y coordinate each point
    
    vector<double> cs;          ///< cosines each point
    vector<double> ss;          ///< sines each point
    vector<double> rs;          ///< radius each point
    vector<double> rfits;       ///< fit radius each point
    
};
