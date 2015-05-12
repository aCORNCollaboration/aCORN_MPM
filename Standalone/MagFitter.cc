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

#include <set>
using std::set;

class FieldDat {
public:
    FieldDat(double mul = 1.0): premul(mul) { }
    
    double max_mul=1.4; ///< maximum multiplier
    double dflt_mul=1;  ///< default multiplier
    double premul;      ///< pre-multiplication of input file for plotting
    TGraph B;           ///< B_z [Gauss] vs. z [m]
    string name;
};

class GordonFieldDat: public FieldDat {
public:
    GordonFieldDat(const string& fname, double mul = 1.0): FieldDat(mul) {
        name = split(split(fname,"/").back(),".")[0];
        string filedat;
        std::ifstream inf;
        inf.open(fname.c_str());
        
        inf >> filedat;
        double sz = 0;
        double sb = 0;
        int igrp = 0;
        size_t nrcols = 6;
        while(inf.good()) {
            auto vdat = sToDoubles(filedat);
            if(vdat.size() < 3*nrcols) continue;
            sz += vdat[0]/100.;
            for(size_t r = 0; r < nrcols; r++) sb += premul*vdat[2+3*r]*10000./nrcols;
            igrp++;
            if(igrp == 1) {
                B.SetPoint(B.GetN(), sz/igrp, sb/igrp);
                sz = sb = igrp = 0;
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
        
        TGraphSmooth gs;
        //B = *(gs.SmoothSuper(&B));
        //B = *(gs.Approx(&B,"linear",nout,xout.data()));
        B = *(gs.SmoothKern(&B, "normal", 0.01, nout, xout.data()));
    
        printf("Loaded %zu field points from '%s'.\n", B.GetN(), fname.c_str());
    }
    
};

struct fitTarget {
    size_t index;       ///< point index
    double B;           ///< target B
    double w;           ///< error weight
};

class FieldShaper {
public:
    FieldShaper() { }
    
    double Bsum(size_t idx) {
        double B = 0;
        for(size_t i=0; i < profiles.size(); i++) {
            if(!zoff[i]) B += profiles[i]->B.GetY()[idx]*currents[i];
            else B += profiles[i]->B.Eval(profiles[i]->B.GetX()[idx] - zoff[i])*currents[i];
        }
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
            double e = it->B - Bsum(it->index);
            assert(e==e);
            sw += it->w;
            serr += e*e*it->w;
        }
        serr /= sw;
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
        ROOT::Math::Functor f(this, &FieldShaper::minFunc, series.size() + moveable.size());
        
        // Choose method upon creation between:
        // kConjugateFR, kConjugatePR, kVectorBFGS,
        // kVectorBFGS2, kSteepestDescent
        // myMin(ROOT::Math::kConjugatePR)
        
        myMin = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
        assert(myMin);
        
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
            myMin->SetLimitedVariable(i, "z_"+to_str(i), 0., 0.001, -0.04, 0.04);
            i++;
        }
        
        printf("\nMinimizing %zu target points over %zu field profiles...\n\n", target.size(), profiles.size());
        
        //fixMovable(true);
        myMin->Minimize();
        //display_parameters(myMin->X());
        
        /*
        fixCurrents(true);
        fixMovable(false);
        myMin->Minimize();
        myMin->Minimize();
        display_parameters(myMin->X());
        
        fixCurrents(false);
        fixMovable(true);
        myMin->Minimize();
        display_parameters(myMin->X());
        */
        
        printf("\nDone!\n\n");
        
        const double* x = myMin->X();        
        display_parameters(x);
        minFunc(x);
        makePlot();
    }
    
    void make_target() {
        assert(profiles.size());
        FieldDat* P = profiles[0];
        
        for(int i=0; i<P->B.GetN(); i++) {
            double z = P->B.GetX()[i];
            fitTarget t;
            t.index = i;
            if(0.25 < z && z < 2.52) {
                t.B = Btarg;
                t.w = (z < 1.1)? 10*(z-0.25)/(1.1-0.25) : (z < 2.2)? 50 : 200;
            } else if(2.70 < z && z < 3.00) {
                    double l = (z-2.7)/(3-2.7);
                    t.B = 220*(1-l)+50*l;
                    t.w = (z < 2.8)? 0.1 : 1e-6;
            } else continue;
            target.push_back(t);
        }
    }
    
    void makePlot() {
        int npts = profiles[0]->B.GetN();
        TGraph gB(npts);
        for(int i = 0; i < npts; i++)
            gB.SetPoint(i, profiles[0]->B.GetX()[i], Bsum(i));
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
    set< set<int> > series;             ///< in-series current groups
    set< set<int> > moveable;           ///< movable coils list, grouped by those moving together
    ROOT::Math::Minimizer* myMin;
    double Btarg = 368.;
    
protected:
    int ncalls = 0;
    vector<double> zoff;                ///< z offsets for each coil
    vector<double> currents;            ///< current (multipliers) for each coil
};


set<int> setone(int i) {
    set<int> s;
    s.insert(i);
    return s;
}

int main(int, char**) {
    
    FieldShaper FS;
    string basedir="/home/mpmendenhall/Documents/aCORN/Reference/CalculatedFieldmaps/";
    FS.profiles.push_back(new GordonFieldDat(basedir+"/Bmap1_20.csv")); // combined 20 coils
    
    FS.EdProf = new GordonFieldDat(basedir+"Bmap_EdSoln_noTrims_I24=14.5.csv");
    
    set<int> B1_21 = setone(FS.profiles.size()-1);
    set<int> B22_23;
    
    for(int i=1; i<=25; i++) {
        set<int> moveset;
        string fname;
        // main
        if(21 <= i && i < 24) {
            string fname = basedir+"/Bmap"+to_str(i);
            //if(22 <= i && i <= 24) fname += "Ed";
            FS.profiles.push_back(new GordonFieldDat(fname+".csv"));
            FS.profiles.back()->B.SetLineColor(2+(i%7));
            moveset.insert(FS.profiles.size()-1);
            
            if(i != 21 && i != 22 && i != 23) FS.series.insert(setone(FS.profiles.size()-1));
            else {
                if(i==21) B1_21.insert(FS.profiles.size()-1);
                else B22_23.insert(FS.profiles.size()-1);
            }
           
        }
        
        // trims
        //fname = basedir+"/Bmap_AxialTrim"+to_str(i);
        fname = basedir+"/Bmap_Trim"+to_str(i);
        //if(22 <= i && i <= 24) fname += "Ed";
        FS.profiles.push_back(new GordonFieldDat(fname+".cvs", 1.));
        moveset.insert(FS.profiles.size()-1);
        FS.profiles.back()->max_mul = 3;
        FS.profiles.back()->dflt_mul = (i==24? 0:0);
        FS.profiles.back()->B.SetLineColor(2+(i%7));
        
        FS.series.insert(setone(FS.profiles.size()-1));
        
        if(i > 22 && i<24) FS.moveable.insert(moveset);
    }
    
    FS.series.insert(B22_23);
    FS.series.insert(B1_21);
    
    FS.make_target();
    FS.optimize_field();
    
    return 0;
}
