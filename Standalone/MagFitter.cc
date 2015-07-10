// We have
// Controlable and stable
// 2ea 8V, 800A supplies
// 1ea 25V 10A that I know works
// 3ea 20V 10A BOPs from Indiana that are flakey, but probably work??
// Several channels 12V 5A
// Many channels 13V 3A
// 
// Cheap and sleazy, not electronically controlable
// 1ea 30V 5A
// 1ea 18V 5A
// 
// 
// Coil Resistance:
// Main Coil ~ 0.18ohm  (approx 44A in one main coil due to voltage limit)
// Trim Coil ~ 3.1 ohm (approx 16 V at 5A)
// 
// We can try to put more than 5A in a trim coil if it is necessary.  5A feels conservative to me, but on the last round that is where we felt we needed to stop for some reason.
//   - Gordon

// make MagFitter; ./MagFitter

#include "PathUtils.hh"
#include "StringManip.hh"
#include "GraphicsUtils.hh"
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphSmooth.h>
#include <TPad.h>
#include <TMath.h>
#include <cassert>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <fstream>
#include <iostream>
#include <algorithm>

class FieldDat {
public:
    FieldDat(double mul = 1.0): premul(mul) { }
    
    double max_mul=1.4; ///< maximum multiplier
    double dflt_mul=1;  ///< default multiplier
    double premul;      ///< pre-multiplication of input file for plotting
    TGraph B;           ///< B_z [Gauss] vs. z [m]
    double I0 = 1.0;    ///< current used to produce profile
    string name;
};

class GordonFieldDat: public FieldDat {
public:
    GordonFieldDat(const string& fname, double mul = 1.0, bool doSmooth = true): FieldDat(mul) {
        name = split(split(fname,"/").back(),".")[0];
        string filedat;
        std::ifstream inf;
        inf.open(fname.c_str());
        if(!inf.good()) {
            printf("Unreadable file: '%s'\n", fname.c_str());
            assert(false);
            return;
        }
        
        inf >> filedat;
        double sz = 0;
        double sb = 0;
        int igrp = 0;
        size_t nrcols = 6;
        vector<double> bsamps;
        while(inf.good()) {
            auto vdat = sToDoubles(filedat);
            if(vdat.size() < 3*nrcols) continue;
            sz += vdat[0]/100.;
            for(size_t r = 0; r < nrcols; r++) {
                bsamps.push_back(premul*vdat[2+3*r]*10000.);
                sb += bsamps.back()/nrcols;
            }
            igrp++;
            if(igrp == 1) {
                B.SetPoint(B.GetN(), sz/igrp, sb/igrp);
                //std::sort(bsamps.begin(), bsamps.end());
                //B.SetPoint(B.GetN(), sz/igrp, bsamps[bsamps.size()/2]);
                sz = sb = igrp = 0;
                bsamps.clear();
            }
            inf >> filedat;
        }
        if(igrp) B.SetPoint(B.GetN(), sz/igrp, sb/igrp);
        
        int nout = 200;
        vector<double> xout;
        double x0 = B.GetX()[0];
        double x1 = B.GetX()[B.GetN()-1];
        for(int i=0; i<nout; i++)
            xout.push_back(x0 + (x1-x0)*i/float(nout-1));
        
        if(doSmooth) {
            TGraphSmooth gs;
            //B = *(gs.SmoothSuper(&B));
            //B = *(gs.Approx(&B,"linear",nout,xout.data()));
            B = *(gs.SmoothKern(&B, "normal", 0.002, nout, xout.data()));
        }
        
        printf("Loaded %zu field points from '%s'.\n", B.GetN(), fname.c_str());
    }
};

struct fitTarget {
    double z;           ///< point location
    double B;           ///< target B
    double w;           ///< error weight
};

class FieldShaper {
public:
    FieldShaper() { }
    
    double Bsum(double z) {
        double B = 0;
        for(size_t i=0; i < profiles.size(); i++) B += profiles[i]->B.Eval(z - zoff[i])*currents[i];
        return B;
    }
    
    double minFunc(const double* params) {
        double serr = 0;
        double sw = 0;
        size_t i = 0;
        for(auto it = series.begin(); it != series.end(); it++) {
            for(auto it2 = it->begin(); it2 != it->end(); it2++)
                currents[*it2] = params[i];
            i++;
        }
        for(auto it = moveable.begin(); it != moveable.end(); it++) {
            for(auto it2 = it->begin(); it2 != it->end(); it2++)
                zoff[*it2] = params[i];
            i++;
        }
        for(auto it = target.begin(); it != target.end(); it++) {
            double e = it->B - Bsum(it->z);
            assert(e==e);
            sw += it->w;
            serr += e*e*it->w;
        }
        serr /= sw;
        
        // diminish shim currents
        for(size_t i=0; i<profiles.size(); i++) {
            double e = currents[i] - profiles[i]->dflt_mul;
            if(!profiles[i]->dflt_mul) serr += 0.001*e*e;
        }
        
        if(!((ncalls++)%20)) printf("Step serr = %g\n", serr);
        return serr;
    }
    
    void display_parameters(const double* x) const {
        size_t i = 0;
        for(auto it = series.begin(); it != series.end(); it++) {
            printf("Coils (");
            //for(auto it2 = it->begin(); it2 != it->end(); it2++) printf(" %i", *it2);
            for(auto it2 = it->begin(); it2 != it->end(); it2++) printf(" %s", profiles[*it2]->name.c_str());
            printf(" )\tcurrent multiplier %g\n", x[i++]);
        }
        for(auto it = moveable.begin(); it != moveable.end(); it++) {
            printf("Moved group (");
            //for(auto it2 = it->begin(); it2 != it->end(); it2++) printf(" %i", *it2);
            for(auto it2 = it->begin(); it2 != it->end(); it2++) printf(" %s", profiles[*it2]->name.c_str());
            printf(" )\t%.02f cm\n", 100*x[i++]);
        }
        
        printf("\n-----------------------\n");
        for(size_t i=0; i<profiles.size(); i++)
            printf("%s:\tI = %+.03f A\n", profiles[i]->name.c_str(), currents[i]*profiles[i]->I0);
        printf("-----------------------\n\n");
    }
    
    void fixCurrents(bool b) {
        for(size_t i=0; i<profiles.size(); i++) {
            if(b)  myMin->FixVariable(i);
            else myMin->ReleaseVariable(i);
        }
    }
    
    void fixMovable(bool b) {
        size_t i = profiles.size();
        for(auto it = moveable.begin(); it != moveable.end(); it++) {
            if(b)  myMin->FixVariable(i++);
            else myMin->ReleaseVariable(i++);
        }
    }
    
    void optimize_field() {
        
        myMin = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
        assert(myMin);
        ROOT::Math::Functor f(this, &FieldShaper::minFunc, series.size() + moveable.size());
        
        myMin->SetMaxFunctionCalls(10000);  // for Minuit
        myMin->SetMaxIterations(10000);     // for GSL
        myMin->SetTolerance(0.001);
        myMin->SetFunction(f);
        
        zoff.resize(profiles.size());
        currents.resize(profiles.size());
        
        size_t i = 0;
        for(auto it = series.begin(); it != series.end(); it++) {
            int j = *(it->begin());
            myMin->SetLimitedVariable(i, "I_"+to_str(i), profiles[j]->dflt_mul, 0.02, -profiles[j]->max_mul, profiles[j]->max_mul);
            i++;
        }
        for(auto it = moveable.begin(); it != moveable.end(); it++) {
            myMin->SetLimitedVariable(i, "z_"+to_str(i), 0., 0.001, -0.06, 0.06);
            i++;
        }
        
        printf("\nMinimizing %zu target points over %zu field profiles...\n\n", target.size(), profiles.size());
        myMin->Minimize();
        printf("\nDone!\n\n");
        
        const double* x = myMin->X();
        minFunc(x);
        display_parameters(x);
        makePlot();
        //errplot((double*)x);
    }
    
    void make_target() {
        assert(profiles.size());
        FieldDat* P = profiles[0];
        
        for(int i=0; i<P->B.GetN(); i++) {
            double z = P->B.GetX()[i];
            fitTarget t;
            t.z = z;
            if(0.25 < z && z < 2.52) {  // collimator = 2.52
                t.B = Btarg;
                t.w = (z < 1.1)? 10*(z-0.25)/(1.1-0.25) : (z < 2.2)? 50 : 200;
            } else if(2.70 < z && z < 3.00) {
                    double l = (z-2.7)/(3-2.7);
                    t.B = 220*(1-l)+50*l;
                    t.w = (z < 2.8)? 0.05 : 1e-6;
            } else continue;
            target.push_back(t);
        }
    }
    
    TGraph makeGraph() {
        int npts = profiles[0]->B.GetN();
        TGraph g(npts);
        for(int i = 0; i < npts; i++) {
            double z = profiles[0]->B.GetX()[i];
            g.SetPoint(i, z, Bsum(z));
        }
        return g;
    }
    
    void errplot(double* x) {
        TGraph g0 = makeGraph();
        x[0] *= 1.001;
        minFunc(x);
        TGraph gp = makeGraph();
        x[0] *= 0.998;
        minFunc(x);
        TGraph gm = makeGraph();
        
        gm.SetLineColor(4);
        gp.SetLineColor(2);
        
        g0.Draw("AC");
        g0.SetTitle("aCORN calculated field profiles");
        g0.GetXaxis()->SetTitle("z position [m]");
        g0.GetXaxis()->SetRangeUser(0,3);
        g0.GetYaxis()->SetTitle("B_{z} [Gauss]");
        g0.SetMinimum(Btarg-5);
        g0.SetMaximum(Btarg+5);
        
        gm.Draw("C");
        gp.Draw("C");
        
        drawHLine(Btarg, gPad, 2, 2);
        drawVLine(2.5135, gPad, 4, 2);
        drawVLine(2.76, gPad, 4, 2);
        gPad->Print("/home/mpmendenhall/Desktop/fieldprof_flux.pdf");
    }
    
    void makePlot() {
        TGraph gB = makeGraph();
        gB.SetMinimum(0);
        gB.Draw("ACP");
        gB.SetTitle("aCORN calculated field profiles");
        gB.GetXaxis()->SetTitle("z position [m]");
        gB.GetXaxis()->SetRangeUser(0,3);
        gB.GetYaxis()->SetTitle("B_{z} [Gauss]");
        //for(auto it = profiles.begin(); it != profiles.end(); it++) (*it)->B.Draw("C");
        if(EdProf) {
            EdProf->B.SetLineColor(2);
            EdProf->B.Draw("C");
        }
        drawHLine(Btarg, gPad, 2, 2);
        drawVLine(2.5135, gPad, 4, 2);
        drawVLine(2.76, gPad, 4, 2);
        
        gPad->Print("/home/mpmendenhall/Desktop/fieldprof.pdf");
        
        gB.SetLineWidth(2);
        gB.Draw("ACP");
        if(EdProf) { EdProf->B.SetLineWidth(2); EdProf->B.Draw("C"); }
        gB.SetMinimum(Btarg-5);
        gB.SetMaximum(Btarg+5);
        gB.GetXaxis()->SetRangeUser(0,3);
        drawHLine(Btarg, gPad, 2, 2);
        drawVLine(2.5135, gPad, 4, 2);
        drawVLine(2.76, gPad, 4, 2);
        gPad->Print("/home/mpmendenhall/Desktop/fieldprof_detail.pdf");
    }
        
    vector<FieldDat*> profiles;         ///< profiles for each coil
    FieldDat* EdProf = NULL;            ///< reference "Ed" profile
    vector<fitTarget> target;           ///< target field profile
    vector< vector<int> > series;       ///< in-series current groups
    vector< vector<int> > moveable;     ///< movable coils list, grouped by those moving together
    ROOT::Math::Minimizer* myMin;
    double Btarg = 368.;
    
protected:
    int ncalls = 0;
    vector<double> zoff;                ///< z offsets for each coil
    vector<double> currents;            ///< current (multipliers) for each coil
};


vector<int> vecone(int i) {
    vector<int> s;
    s.push_back(i);
    return s;
}

int main(int, char**) {
    
    FieldShaper FS;
    string basedir="/home/mpmendenhall/Documents/aCORN/Reference/CalculatedFieldmaps/";
    FS.profiles.push_back(new GordonFieldDat(basedir+"/Bmap1_20.csv")); // combined 20 coils
    FS.profiles.back()->I0 = 29.2;
    
    FS.EdProf = new GordonFieldDat(basedir+"Bmap_EdSoln_noTrims_I24=14.5.csv", 1.0, false);
    
    vector<int> Bmain = vecone(FS.profiles.size()-1);
    vector<int> B22_23;
    
    bool EdSpacing = false;
    
    for(int i=1; i<=25; i++) {
        vector<int> moveset;
        string fname;
        // main
        if(21 <= i && i <= 23+EdSpacing) {
            string fname = basedir+"/Bmap"+to_str(i);
            if(EdSpacing && 22 <= i && i <= 24) fname += "Ed";
            else if(i==23) fname = basedir+"/Bmap_Main23-6cm";
            FS.profiles.push_back(new GordonFieldDat(fname+".csv"));
            FS.profiles.back()->B.SetLineColor(2+(i%7));
            FS.profiles.back()->I0 = 29.2;
            moveset.push_back(FS.profiles.size()-1);
            
            if(i <= 21) Bmain.push_back(FS.profiles.size()-1);
            else B22_23.push_back(FS.profiles.size()-1);
            //else FS.series.push_back(vecone(FS.profiles.size()-1));
           
        }
        
        // trims
        if(i==12 && !EdSpacing) continue; // turn off middle trim to force reduced trim currents
        
        fname = basedir+"/Bmap_Trim"+to_str(i);
        if(EdSpacing && 21 <= i && i <= 25) fname += "Ed";
        else if(i==23) fname = basedir+"/Bmap_Trim23-6cm";
        FS.profiles.push_back(new GordonFieldDat(fname+".csv"));
        moveset.push_back(FS.profiles.size()-1);
        FS.profiles.back()->max_mul = 3;
        FS.profiles.back()->dflt_mul = (i==24? 0:0);
        FS.profiles.back()->B.SetLineColor(2+(i%7));
        
        FS.series.push_back(vecone(FS.profiles.size()-1));
        
        //if(i == 23 || i == 22) FS.moveable.push_back(moveset);
    }
    
    if(B22_23.size() >= 2) FS.series.push_back(B22_23);
    FS.series.push_back(Bmain);
    
    FS.make_target();
    FS.optimize_field();
    
    return 0;
}
