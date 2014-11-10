#include "NuclEvtGen.hh"
#include "Collimator.hh"
#include "SourceCalPlugin.hh"
#include "PathUtils.hh"

#include <cmath>
#include <map>

#include <Math/QuasiRandom.h>
#include <TStyle.h>
#include <TMath.h>
#include <TGraph.h>

using namespace ROOT::Math;
using std::vector;
using std::pair;
using std::map;

/// simple calculations for electron survival through collimator
class ElectronSurvival {
public:
    /// constructor
    ElectronSurvival() { }
    
    double B0 = 364;    ///< B field [Gauss]
    double l = 2.75;    ///< electron collimator radius [cm]
    double d = 0;       ///< center offset
    double m_e = 511.;  ///< electron mass
    
    /// survival fraction for given kinetic energy
    double survival_frac(double KE) {
        if(precalc.count(KE)) return precalc[KE];
        
        if(d > l) return 0.;
        double a0 = 3.34 * sqrt(KE*KE+2*m_e*KE)/B0; // max Larmor radius
        if(d + 2*a0 <= l) return 1.;
        
        // cos phi below which no trajectories survive
        double c0 = (l+d < 2*a0? sqrt(1-pow((l+d)/(2.*a0),2)) : 0);
        // cos phi above which all trajectories survive
        double c1 = sqrt(1-pow((l-d)/(2.*a0),2));
   
        double f = nintegrate_survival(a0, c0, c1) + 1-c1;
        precalc[KE] = f;
        printf("E=%g: survival = %g\n", KE, f);
        return f;
    }
    
    map<double,double> precalc;
    
protected:
    
    /// fraction surviving as a function of cos(phi)
    double _survival_frac(double a0, double c) const {
        double a = a0*sqrt(1-c*c);
        if(d >= l or 2*a >= d+l) return 0.;
        if(d + 2*a <= l) return 1.;
        if(d == 0) return 2*a < l;
        double x = (l*l - d*d - 2*l*a)/(2*d*a);
        return x<=-1? 0 : (x >= 1? 1 : 1 - acos(x)/M_PI);
    }
    
    /// numerically integrate survival fraction (Simpson's Rule)
    double nintegrate_survival(double a0, double x1, double x2, int npts=100) const {
        double s = 0;
        for(int i=0; i<=npts; i++) {
            double x = x1 + i*(x2-x1)/npts;
            double w = i%2? 4 : 2;
            if(i==0 || i==npts) w = 1;
            s += _survival_frac(a0, x)*w;
        }
        return s*(x2-x1)/(3.*npts);
    }
};

/// Scintillator quenching estimator        
class QuenchCalculator {
public:
    /// Constructor
    QuenchCalculator(size_t npts = 1000): gQ(npts) {
        for(size_t i = 0; i<npts; i++) {
            double Edep = i*2000./(npts-1);
            gQ.SetPoint(i,Edep,Equench(Edep));
        }
    }
    
    double n0 = 0.5426;
    double rho = 1.0;
    double I = 62.5e-6;
    double kB = 0.015;
    double Emirror = 3.0;     ///< electrostatic mirror energy boost [keV]
    
    /// Bethe formula for energy deposition
    double bethe(double KE) const {
        if(KE <= 0) return 0;
        double b = sqrt(1.-1./pow(KE/511.+1,2));
        double bth = n0 * rho * 0.307/(b*b)*(log(1.022*b*b/(1-b*b)/I)-b*b);
        return bth>0? bth : 0;
    }
    
    double EQfromEdep(double Edep) const { return gQ.Eval(Edep+Emirror)/0.945; }
    
protected:
    TGraph gQ;  ///< quenching curve; x=E_dep, y=E_Q
    
    /// calculate quenched energy
    double Equench(double Edep, int npts=100) const {
        double s = 0;
        for(int i=0; i<=npts; i++) {
            double x = i*Edep/npts;
            double w = i%2? 4 : 2;
            if(i==0 || i==npts) w = 1;
            s += w/(1+kB*bethe(x));
        }
        return s*Edep/(3.*npts);
    }
};

/// Fill histogram preserving *average* value interpolated between bins
void fill_interp(TH1* h, double x, double w = 1.0) {
    TAxis* Ax = h->GetXaxis();
    int b0 = Ax->FindBin(x);
    double c0 = Ax->GetBinCenter(b0);
    int b1 = x > c0? b0+1 : b0-1;
    double c1 = Ax->GetBinCenter(b1);
    double a = (c1-x)/(c1-c0);
    h->Fill(c0, a*w);
    h->Fill(c1, (1-a)*w);
}

/// Fill all probabilistic combinations of energies
void fill_e_combos(double E0, double p0, pair<double,double>* ee, size_t n, TH1* h) {
    if(!n) { if(E0) fill_interp(h,E0,p0); }
    else {
        fill_e_combos(E0+ee[0].first, p0*ee[0].second, ee+1, n-1, h);
        fill_e_combos(E0, p0*(1. - ee[0].second), ee+1, n-1, h);
    }
}


/// Poisson-ize energy histogram; preserving total counts
void poisson_smear(const TH1& hIn, TH1& hOut, double PEperKeV, QuenchCalculator* QC = NULL) {
    for(int i=1; i<=hIn.GetNbinsX(); i++) {
        double c0 = hIn.GetBinContent(i);
        if(!c0) continue;
        double E0 = hIn.GetBinCenter(i);
        double n0 = (QC?QC->EQfromEdep(E0):E0)*PEperKeV;
        double nrm = 0;
        for(int j=1; j<=hOut.GetNbinsX(); j++) nrm += TMath::Poisson(hOut.GetBinCenter(j)*PEperKeV, n0);
        for(int j=1; j<=hOut.GetNbinsX(); j++) {
            double E1 = hOut.GetBinCenter(j);
            hOut.Fill(E1, c0 * TMath::Poisson(E1*PEperKeV, n0)/nrm);
        }
    }
}

// generate calibration spectra
void gen_calib_spectra(const std::string& sname, TH1*& hInitEnergy, TH1*& hEnergy) {
    
    NucDecayLibrary NDL(getEnvSafe("ACORN_AUX")+"/NuclearDecays/");
    NucDecaySystem& srcgen = NDL.getGenerator(sname);
    
    ElectronSurvival ES;
    
    QuasiRandomSobol QR(srcgen.getNDF());
    double u[128];
    
    hInitEnergy = new TH1F("hInitEnergy", (sname + " electrons").c_str(), 1000, 0, 2000);
    hInitEnergy->GetXaxis()->SetTitle("Energy [keV]");
    hInitEnergy->GetYaxis()->SetTitle("Counts / keV / 1e6 decays");
    hInitEnergy->GetYaxis()->SetTitleOffset(1.4);
    
    hEnergy = new TH1F("hEnergy", (sname + " electrons").c_str(), 2000, 0, 2000);
    hEnergy->GetXaxis()->SetTitle("Energy [keV]");
    hEnergy->GetYaxis()->SetTitleOffset(1.4);
    
    const size_t npts = 1000000;
    for(size_t i=0; i<npts; i++) {
        QR.Next(u);
        vector<NucDecayEvent> v;
        srcgen.genDecayChain(v, u);
        
        vector< pair<double,double> > eevts;
        for(auto it = v.begin(); it != v.end(); it++) {
            if(it->d == D_ELECTRON) {
                hInitEnergy->Fill(it->E);
                eevts.push_back(pair<double,double>(it->E, 0.5*ES.survival_frac(it->E)));
            }
        }
        fill_e_combos(0., 1., &eevts[0], eevts.size(), hEnergy);
    }
    
    hInitEnergy->Scale(1e6/npts);
    hEnergy->Scale(1e6/npts);
}

// load calibration spectra from file
void load_calib_spectra(const std::string& sname, TH1*& hInitEnergy, TH1*& hEnergy) {
    
    TFile felectron((getEnvSafe("ACORN_MCOUT") + "/aCORN_" + sname + "_e/Plots/aCORN_Scatter.root").c_str(),"READ");
    hInitEnergy = (TH1*)felectron.Get("hEnergy0");      // counts / 1e6 decays
    hEnergy = (TH1*)felectron.Get("hEnergy1");          // counts / 1e6 decays
    assert(hEnergy && hInitEnergy);
    
    TFile fgamma((getEnvSafe("ACORN_MCOUT") + "/aCORN_" + sname + "_gamma/Plots/aCORN_Scatter.root").c_str(),"READ");
    TH1* hInitEnergyG = (TH1*)fgamma.Get("hEnergy0");   // counts / 1e6 decays
    TH1* hEnergyG = (TH1*)fgamma.Get("hEnergy1");       // counts / 1e6 decays
    if(hInitEnergyG && hEnergyG) {
        double gammaScale = 3.5;
        hInitEnergy->Add(hInitEnergyG, gammaScale);
        hEnergy->Add(hEnergyG, gammaScale);
    } else {
        printf("Gamma constributions not found for '%s'\n", sname.c_str());
    }
}

// load data spectrum
TH1* loadData() {
    TFile f((getEnvSafe("ACORN_SUMMARY")+"SourceCal/SourceCal.root").c_str(),"READ");
    TH1* hdat = (TH1*)f.Get("hEnergy_Rate");
    assert(hdat);
    return hdat;
}

int main(int, char**) {
    
    const std::string& sname = "Bi207";
    
    OutputManager OM("Simulated", getEnvSafe("ACORN_SUMMARY")+"/Sim_"+sname+"/");

    TH1* hInitEnergy;
    TH1* hEnergy;
    //gen_calib_spectra(sname, hInitEnergy, hEnergy);
    load_calib_spectra(sname, hInitEnergy, hEnergy);
    
    TH1F hSmeared("hSmeared", (sname + " PE-smeared energy").c_str(), 500, 0, 1500);
    hSmeared.GetXaxis()->SetTitle("Quenched Energy [keV]");
    hSmeared.GetYaxis()->SetTitle("Counts / 1e6 decays");
    hSmeared.GetYaxis()->SetTitleOffset(1.4);

    
    gStyle->SetOptStat("");
    OM.defaultCanvas->SetLogy(true);
    hInitEnergy->Draw("HIST");
    hEnergy->Draw("HIST SAME");
    OM.printCanvas("Energy");
    OM.defaultCanvas->SetLogy(false);
    
    QuenchCalculator QC;
    poisson_smear(*hEnergy, hSmeared, 0.22, &QC);
    hSmeared.SetMaximum(3);
    
    TF1 fGaus("fGaus","gaus",400,550);
    hSmeared.Fit(&fGaus,"R");
    
    hSmeared.Draw("HIST");
    OM.printCanvas("EnergyPE");
    
    TH1* hdat = loadData();
    hdat->Draw();
    hSmeared.Scale(26.);
    hSmeared.Draw("HIST SAME");
    OM.printCanvas("DataComparison");

    
    return 0;
}

/*
gSystem->Load("libEventLib.so")
TFile f("scint_Bi207_gammas.root"); TTree* T = f.Get("PG4")
TH1F h("hGammas","gamma energy deposition",2000,0,2000)
T.Draw("1000*EIoni >> hGammas")
h.Scale(1./557.33)
h.GetXaxis()->SetTitle("deposited energy [keV]")
h.GetYaxis()->SetTitle("counts / keV / 10000 decays")
TFile fout("Bi207_Gammas.root","RECREATE")
fout.cd()
h.Write()
*/
