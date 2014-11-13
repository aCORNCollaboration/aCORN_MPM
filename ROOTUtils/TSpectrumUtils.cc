#include "TSpectrumUtils.hh"
#include "PathUtils.hh"
#include "GraphicsUtils.hh"

bool yrevsort(std::pair<float,float> a, std::pair<float,float> b) { return (b.second < a.second); }

vector<TSpectrumPeak> tspectrumSearch(TH1* hin, float sigma, float thresh) {
    TSpectrum* TS = new TSpectrum();
    int npks = TS->Search(hin,sigma,"",thresh);
    if(npks<0)
        npks = 0;
    
    vector<TSpectrumPeak> vx;
    return vx;
}

#ifdef TSPECTRUM_USES_DOUBLE
typedef Double_t TSpectrum_Data_t;
#else
typedef Float_t TSpectrum_Data_t;
#endif

vector<float> tspectrumSearch(TH1* hin, TH1** hout, float sigma, float thresh) {
    
    // temporary data arrays
    unsigned int nbins = hin->GetNbinsX();
    TSpectrum_Data_t* datin = new TSpectrum_Data_t[nbins];
    TSpectrum_Data_t* datout = new TSpectrum_Data_t[nbins];
    for(unsigned int i=0; i<nbins; i++)
        datin[i] = hin->GetBinContent(i+1);
    
    // TSpectrum fitting
    TSpectrum* TS = new TSpectrum();
    int npks = TS->SearchHighRes(datin,datout,nbins,sigma,thresh,true,1,false,1);
    if(npks<0)
        npks = 0;
    
    // output peak positions vector
    vector<float> vx;
    vector< std::pair<float,float> > fndpks;
    for(int i=0; i<npks; i++) {
        float binpos = TS->GetPositionX()[i];
        fndpks.push_back(std::make_pair(binterpolate(hin->GetXaxis(),binpos),datout[int(binpos+0.5)]));
    }
    std::sort(fndpks.begin(),fndpks.end(),yrevsort);
    for(int i=0; i<npks; i++)
        vx.push_back(fndpks[i].first);
    
    // output histogram
    if(hout) {
        *hout = (TH1*)hin->Clone("hout");
        for(unsigned int i=0; i<nbins; i++)
            (*hout)->SetBinContent(i+1,datout[i]);
        (*hout)->SetBinContent(0,0);
        (*hout)->SetBinContent(nbins+1,0);
        (*hout)->SetLineColor(2);
    }
    
    // cleanup
    delete[](datin);
    delete[](datout);
    delete(TS);
    
    return vx;
}

vector<SpectrumPeak> tspectrumPrefit(TH1* indat, float searchsigma, const vector<SpectrumPeak>& expectedPeaks,
                                     TH1*& htout, float pkMin, float pkMax) {
    
    // use TSpectrum peak fitting to find initial guesses for peak locations
    float binsigma = searchsigma/indat->GetBinWidth(1);
    if(binsigma < 2.) {
        searchsigma *= 2./binsigma;
        binsigma = 2.;
    }
    vector<float> tpks0 = tspectrumSearch(indat,&htout,binsigma,10.0);
    
    // throw out out-of-range peaks
    vector<float> tpks;
    for(vector<float>::iterator it = tpks0.begin(); it != tpks0.end(); it++)
        if(pkMin <= *it && *it <= pkMax)
            tpks.push_back(*it);
        
        // limit to number of found peaks
        vector<SpectrumPeak> foundPeaks = expectedPeaks;
    if(tpks.size() < expectedPeaks.size()) {
        foundPeaks.clear();
        return foundPeaks;
    }
    while(foundPeaks.size()>tpks.size())
        foundPeaks.pop_back();
    unsigned int npks = foundPeaks.size();
    
    // pre-fit the TSpectrum peak-by-peak to estimate MultiGaus parameters
    vector<TF1> tspks;
    for(unsigned int n=0; n<npks; n++) {
        tspks.push_back(TF1("PeakFit","gaus",tpks[n]-searchsigma,tpks[n]+searchsigma));
        tspks.back().SetParameter(0,indat->GetBinContent(indat->FindBin(tpks[n])));
        tspks.back().FixParameter(1,tpks[n]);
        tspks.back().SetParameter(2,searchsigma);
        htout->Fit(&tspks.back(),"QR+");
        foundPeaks[n].fromGaussian(&tspks.back());
    }
    
    return foundPeaks;
}
                                     
MultiGaus multiPeakFitter(TH1* indat, const vector<SpectrumPeak>& expectedPeaks, float nSigma) {
    unsigned int npks = expectedPeaks.size();
    MultiGaus mg(npks,"SpectrumFitter",nSigma);
    for(unsigned int n=0; n<npks; n++) {
        mg.setParameter(3*n+0, expectedPeaks[n].h.x);
        mg.setParameter(3*n+1, expectedPeaks[n].center.x);
        mg.setParameter(3*n+2, expectedPeaks[n].width.x);
    }
    mg.getFitter()->SetLineColor(1); // 6
    mg.fit(indat);
    mg.fit(indat);
    
    return mg;
}

vector<SpectrumPeak> fancyMultiFit(TH1* indat, float searchsigma, const vector<SpectrumPeak>& expectedPeaks,
                                   bool bgsubtract, const std::string& drawName, float nSigma, float pkMin, float pkMax) {
    
    // use TSpectrum peak fitting to find initial guesses for peak locations
    TH1* htout;
    vector<SpectrumPeak> foundPeaks = tspectrumPrefit(indat, searchsigma, expectedPeaks, htout, pkMin, pkMax);
    
    if(bgsubtract) {
        TH1F* bgspec = (TH1F*)TSpectrum().Background(indat);
        indat->Add(bgspec,-1.0);
        delete(bgspec);
    }
    
    Size_t npks = foundPeaks.size();
    if(npks == expectedPeaks.size()) {
        MultiGaus mg = multiPeakFitter(indat, foundPeaks, nSigma);
        
        // extract data
        for(unsigned int p=0; p<npks; p++) {
            foundPeaks[p].h = mg.getPar(3*p+0);
            float_err x0 = foundPeaks[p].energyCenter = mg.getPar(3*p+1);
            float_err w = foundPeaks[p].energyWidth = mg.getPar(3*p+2);
            foundPeaks[p].integral = indat->Integral(indat->FindBin(x0.x-w.x),indat->FindBin(x0.x+w.x))*indat->GetBinWidth(1);
            foundPeaks[p].center = foundPeaks[p].width = 0; // zero out old pre-fit estimates
        }
        
        // draw histograms
        if(drawName.size()) {
            makePath(drawName,true);
            indat->Draw();
            //htout->Draw("Same");
            for(unsigned int p=0; p<npks; p++) {
                drawVLine(foundPeaks[p].energyCenter.x, gPad, 1);
                drawVLine(foundPeaks[p].energyCenter.x-foundPeaks[p].energyWidth.x, gPad, 1, 2);
                drawVLine(foundPeaks[p].energyCenter.x+foundPeaks[p].energyWidth.x, gPad, 1, 2);
            }
            gPad->Print(drawName.c_str());
        }
        
        // remove fit from histogram associated functions
        indat->GetListOfFunctions()->Remove(indat->GetListOfFunctions()->FindObject(mg.getFitter()->GetName()));
        
        delete(htout);
        return foundPeaks; 
    } else {
        // draw histograms
        if(drawName.size()) {
            indat->Draw();
            htout->Draw("Same");
            gPad->Print(drawName.c_str());
        }
        delete(htout);
        return foundPeaks;
    }
}
