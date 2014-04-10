#ifndef HISTOGRAMSEQUENCEFITTER_HH
#define HISTOGRAMSEQUENCEFITTER_HH

#include <cmath>
#include <vector>

#include <TH1.h>
#include <TF1.h>
#include <TF2.h>
#include <TFitResultPtr.h>
#include <TGraph2DErrors.h>

/// Bin-by-bin fits for histograms integrating over a time-evolving quantity
class HistogramSequenceFitter {
public:
	/// constructor
	HistogramSequenceFitter(unsigned int npar) { f = new TF2("f", this, 0, 1, 0, 1, npar, "HistogramSequenceFitter"); }
	/// destructor
	virtual ~HistogramSequenceFitter() {}
	
	/// add fit data point
	void addData(const TH1* h, double t1, double t2);
	/// calculate fit results
	void fit();
	/// generate interpolated histogram for time range
	TH1* interpolate(double t1, double t2);
	
	/// fitting evaluation function
	double operator() (double *x, double *p) { return integ_f(x[1],p)-integ_f(x[0],p); }
	
protected:
	
	/// indefinite integral of variation function (subclass me!)
	virtual double integ_f(double t, const double* params) const = 0;
	/// initial fit parameter estimation
	virtual void estimate_params(const TGraph2DErrors& d) = 0;
	
	TF2* f;	 						//< fitter function operating on this class
	std::vector<const TH1*> hs;		//< input data point histograms
	std::vector<double> T1;			//< start times
	std::vector<double> T2;			//< end times
	
	std::vector<TFitResultPtr> fts;	//< fits for each bin
};

/// Exponentially decaying sequence
class ExponentialSequenceFitter: public HistogramSequenceFitter {
public:
	/// constructor
	ExponentialSequenceFitter(): HistogramSequenceFitter(2), T0(0) { }
	
protected:

	/// indefinite integral of variation function (subclass me!)
	virtual double integ_f(double t, const double* params) const { return params[0] * exp((t-T0)*params[1]); }
	/// initial fit parameter estimation
	virtual void estimate_params(const TGraph2DErrors& d);
	
	static TF1 expFit;	//< fitter for parameter estimation
	double T0;			//< fit reference start time
};

/// Polynomially varying sequence
class PolynomialSequenceFitter: public HistogramSequenceFitter {
public:
	/// constructor
	PolynomialSequenceFitter(unsigned int npar): HistogramSequenceFitter(npar), T0(0) { }
	
protected:

	/// indefinite integral of variation function (subclass me!)
	virtual double integ_f(double t, const double* params) const;
	/// initial fit parameter estimation
	virtual void estimate_params(const TGraph2DErrors& d);
	
	double T0;			//< fit reference start time
};

#endif
