#ifndef TSPECTRUMUTILS_HH
#define TSPECTRUMUTILS_HH 1

#include <vector>
#include <TH1F.h>
#include <TPad.h>
#include <TSpectrum.h>
#include "MultiGaus.hh"
#include "SpectrumPeak.hh"
#include "Types.hh"
#include <float.h>
#include <algorithm>

/// simple struct for peaks identified by TSpectrum
struct TSpectrumPeak {
    float x; ///< peak position
    float y; ///< peak height
};

/// more convenient wrapper for TSpectrumSearch
vector<TSpectrumPeak> tspectrumSearch(TH1* hin, float sigma = 2.0, float thresh = 0.01);

/// perform a TSpectrum peak search on the input histogram for peaks of width sigma bins; report results up to specified threshold below top peak
vector<float> tspectrumSearch(TH1* hin, TH1** hout, float sigma, float thresh);

/// pre-fit using a TSpectrum to estimate peak positions
vector<SpectrumPeak> tspectrumPrefit(TH1* indat, float searchsigma, const vector<SpectrumPeak>& expectedPeaks,
                                     TH1*& htout, float pkMin = -FLT_MAX, float pkMax = FLT_MAX);

/// generate fitter for fitting multi-peak spectum
MultiGaus multiPeakFitter(TH1* indat, const vector<SpectrumPeak>& expectedPeaks, float nSigma = 1.5);

/// locate peaks in multi-peak spectrum based on initial guesses
vector<SpectrumPeak> fancyMultiFit(TH1* indat, float searchsigma, const vector<SpectrumPeak>& expectedPeaks,
                                   bool bgsubtract = false, const std::string& drawName = "", float nSigma = 1.5, float pkMin = -FLT_MAX, float pkMax = FLT_MAX);

#endif
