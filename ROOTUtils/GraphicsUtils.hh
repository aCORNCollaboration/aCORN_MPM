#ifndef GRAPHICSUTILS_HH
#define GRAPHICSUTILS_HH

#include <math.h>
#include <TVirtualPad.h>
#include "Enums.hh"

#include <vector>
#include <TH1.h>
#include <cfloat>

using std::vector;

/// draw several histograms simultaneously; return max histogram height
double drawSimulHistos(vector<TH1*>& hists, const std::string& opt = "", const std::string& newTitle = "DEFAULT");
/// draw a pair of histograms (red and blue by default)
void drawHistoPair(TH1* hRed, TH1* hBlue, const std::string& opt = "", Int_t c1 = 2, Int_t c2 = 4);
/// draw a pair of histograms, with line and diamond markers
void drawDataMCPair(TH1* dat, TH1* mc);

/// draw vertical line marker
void drawVLine(Float_t x, TVirtualPad* C, Int_t color = 4, Int_t style = 1);
/// draw horizontal line
void drawHLine(Float_t y, TVirtualPad* C, Int_t color = 4);
/// draw circle
void drawCircle(float r, Int_t color = 1, Int_t lstyle = 1, float x0=0, float y0=0);
/// draw shaded rectangle marker
void drawExcludedRegion(Float_t x0, Float_t x1, TCanvas* C, Int_t color = 4, Int_t fill = 1001);

/// set up grayscale figures color palette, black-to-white by default, or inverted
void makeGrayscalepalette(bool b2w = true);

/// set up blue-white-red palette
void makeRBpalette();

#endif

