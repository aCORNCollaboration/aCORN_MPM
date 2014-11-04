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
    
    /// Bethe formula for energy deposition
    double bethe(double KE) const {
        if(KE <= 0) return 0;
        double b = sqrt(1.-1./pow(KE/511.+1,2));
        double bth = n0 * rho * 0.307/(b*b)*(log(1.022*b*b/(1-b*b)/I)-b*b);
        return bth>0? bth : 0;
    }
    
    double EQfromEdep(double Edep) const { return gQ.Eval(Edep); }
    
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

/// Poisson-ize energy histogram
void poisson_smear(const TH1& hIn, TH1& hOut, double PEperKeV) {
    for(int i=1; i<=hIn.GetNbinsX(); i++) {
        double c0 = hIn.GetBinContent(i);
        if(!c0) continue;
        double n0 = hIn.GetBinCenter(i)*PEperKeV;
        for(int j=1; j<=hOut.GetNbinsX(); j++) {
            double E1 = hOut.GetBinCenter(j);
            hOut.Fill(E1, c0 * TMath::Poisson(E1*PEperKeV, n0));
        }
    }
}

int main(int, char**) {
    
    const std::string& sname = "Bi207";
    
    OutputManager OM("Simulated", getEnvSafe("ACORN_SUMMARY")+"/Sim_"+sname+"/");
    
    NucDecayLibrary NDL(getEnvSafe("ACORN_AUX")+"/NuclearDecays/");
    NucDecaySystem& srcgen = NDL.getGenerator(sname);
    
    ElectronSurvival ES;
    QuenchCalculator QC;
    
    QuasiRandomSobol QR(srcgen.getNDF());
    double u[128];
    
    TH1F hInitEnergy("hInitEnergy", (sname + " electrons").c_str(), 1000, 0, 2000);
    hInitEnergy.GetXaxis()->SetTitle("Energy [keV]");
    hInitEnergy.GetYaxis()->SetTitle("Counts per 1000 decays");
    hInitEnergy.GetYaxis()->SetTitleOffset(1.4);
    
    TH1F hEnergy("hEnergy", (sname + " electrons").c_str(), 2000, 0, 2000);
    hEnergy.GetXaxis()->SetTitle("Energy [keV]");
    TH1F hSmeared("hSmeared", (sname + " PE-smeared energy").c_str(), 500, 0, 1500);
    hSmeared.GetXaxis()->SetTitle("Quenched Energy [keV]");
    
    const size_t npts = 1000000;
    for(size_t i=0; i<npts; i++) {
        QR.Next(u);
        vector<NucDecayEvent> v;
        srcgen.genDecayChain(v, u);
        
        vector< pair<double,double> > eevts;
        for(auto it = v.begin(); it != v.end(); it++) {
            if(it->d == D_ELECTRON) {
                hInitEnergy.Fill(it->E);
                eevts.push_back(pair<double,double>(QC.EQfromEdep(it->E), 0.5*ES.survival_frac(it->E)));
            }
        }
        fill_e_combos(0., 1., &eevts[0], eevts.size(), &hEnergy);
    }
    
    double ntot = hEnergy.Integral();
    hEnergy.Scale(1./ntot);
    
    gStyle->SetOptStat("");
    OM.defaultCanvas->SetLogy(true);
    hEnergy.SetMinimum(1e-4);
    hEnergy.SetMaximum(1);
    hEnergy.Draw("HIST");
    OM.printCanvas("Energy");
    
    hInitEnergy.Scale(1000./npts);
    hInitEnergy.Draw("HIST");
    OM.printCanvas("InitEnergy");
    
    // lousy gamma tail estimate
    //for(int j=1; j<=hSmeared.GetNbinsX(); j++) {
    //    double E1 = hSmeared.GetBinCenter(j);
    //    hSmeared.SetBinContent(j, 0.1/E1);
    //}
    poisson_smear(hEnergy, hSmeared, 0.26);
    hSmeared.Scale(1000);
    
    TF1 fGaus("fGaus","gaus",400,550);
    hSmeared.Fit(&fGaus,"R");
    
    OM.defaultCanvas->SetLogy(false);
    hSmeared.SetMinimum(0);
    hSmeared.SetMaximum(2.5);
    hSmeared.Draw("HIST");
    OM.printCanvas("EnergyPE");
    
    return 0;
}
