/// \file FieldMapFitter.cc aCORN magnetic field map calculations
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

// make FieldMapFitter; ./FieldMapFitter

#include "PathUtils.hh"
#include "StringManip.hh"
#include "GraphicsUtils.hh"

#include <cassert>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include <TMath.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

/// Transform of coordinates by rotation around bore axis
class RotationTransform {
public:
    /// Constructor
    RotationTransform(double th): theta(th), s(sin(theta)), c(cos(theta)) { }
    
    void apply(const double xin[3], double xout[3]) {
        xout[2] = xin[2];
        xout[0] = c*xin[0] - s*xin[1];
        xout[1] = s*xin[0] + c*xin[1];
    }
    
    double theta, s, c;
};

/// parametrized B Field model
class BModel {
public:
    BModel() { }
    
    virtual double fieldAtAlong(const double x[3], const double v[3], const double* p) const {
        double B[3] = { p[0] + p[3]*x[0] + p[4]*x[1], 
                        p[1] + p[5]*x[0] + p[6]*x[1],
                        p[2] + p[7]*x[0] + p[8]*x[1] };
        return B[0]*v[0] + B[1]*v[1] + B[2]*v[2];
    }
    
    virtual void display(const double* p) const {
        printf("B = < %g\t%g\t%g > Gauss + gradients", p[0], p[1], p[2]);
        for(int i=3; i<9; i++) printf("\t%g",p[i]);
        printf("\n");
    }
    
    virtual int nparams() const { return 9; }
};

class ConstBModel: public BModel {
public:
    ConstBModel() { }
    
    double fieldAtAlong(const double[3], const double v[3], const double* p) const override {
        return p[0]*v[0] + p[1]*v[1] + p[2]*v[2];
    }
    
    void display(const double* p) const override {
        printf("B = < %g\t%g\t%g > Gauss\n", p[0], p[1], p[2]);
    }
    
    int nparams() const override { return 3; }
    
};

/// Fitter class for field map data at fixed z
class MapFitter {
public:
    /// Constructor
    MapFitter() { }
    
    vector<RotationTransform> rots;
    
    double scan_z;      ///< Z position [cm]
    double scan_t;      ///< scan time [s]
    double scan_T;      ///< temperature [deg. C]
    bool offcenter = false;     ///< whether using off-center probes
    
    struct meas_point {
        int p_idx;      ///< probe index
        int t_idx;      ///< transform index
        double val;     ///< measurement value (B field component in Gauss)
    };
    vector<meas_point> meas;
    
    struct probe_coords {
        double x0[3];   ///< un-transformed position
        double v0[3];   ///< un-transformed direction (unit vector)
    };
    vector<probe_coords> probes;
    
    ConstBModel myB;
    //BModel myB;
    int ncalls = 0;
    
    /// minimization criterea
    double minFunc(const double* params) {
        double serr = 0; // sum of error terms
        
        // set probe directions from parameters
        const double* p = params + myB.nparams();
        for(size_t n=0; n < probes.size(); n++) {
            double pnorm = 0;
            for(int i=0; i<3; i++) { 
                probes[n].v0[i] = p[3*n+i]; 
                pnorm += probes[n].v0[i] * probes[n].v0[i];
            }
            pnorm = sqrt(pnorm);
            assert(pnorm);
            for(int i=0; i<3; i++) probes[n].v0[i] /= pnorm;
            serr += (pnorm-1)*(pnorm-1); // normalization constraints
        }
        
        // evaluate fields for each measurement
        double x1[3];
        double v1[3];
        for(auto& m: meas) {
            auto& pb = probes[m.p_idx];
            auto& tx = rots[m.t_idx];
            tx.apply(pb.v0, v1);
            tx.apply(pb.x0, x1);
            double dB = m.val - myB.fieldAtAlong(x1, v1, params);
            serr += dB * dB;
        }
        
        return serr;
    }

    /// set up rotations for equal-angle circular scans
    void setup_circle_rots(int n) {
        rots.clear();
        for(int i=0; i<n; i++) rots.push_back(RotationTransform(i*2*M_PI/n));
    }
    
    /// set up probes
    void set_probes() {
        if(offcenter) {
            probes.resize(2); // z probe is dead.
            for(int i=0; i<2; i++) {
                for(int a=0; a<3; a++) {
                    probes[i].x0[a] = (a==0)? 5 : 0;
                    probes[i].v0[a] = (a==i)? (i?-1:1) : 0;
                }
            }
        } else {
            probes.resize(3);
            for(int i=0; i<3; i++) {
                for(int a=0; a<3; a++) {
                    probes[i].x0[a] = 0;
                    probes[i].v0[a] = (a==i)? 1 : 0;
                }
            }
        }
        display_probes();
    }
    
    /// load measurement data line for n circular points (n=4 or n=12); includes extra wrap-around point
    void load_circle_data(const string& l, size_t n) {
        assert(!(12%n));
        assert(n == rots.size());
        int ngroup = 12/n;
        
        auto v = sToDoubles(l);
        scan_z = v[0];
        scan_t = v[1];
        scan_T = v[2];
        size_t np = probes.size();
        meas.resize(np*(n+1));
        
        printf("Field points:\n");
        for(size_t i=0; i<n+1; i++) {
            printf("<\t");
            for(size_t a=0; a<probes.size(); a++) {
                meas[np*i+a].p_idx = a;
                meas[np*i+a].t_idx = i%n;
                meas[np*i+a].val = v[3 + 4*ngroup*i + a];
                printf("%g\t",  meas[np*i+a].val);
            }
            printf(">\n");
        }
    }
    
    /// display probe (fitted) orientations
    void display_probes() const {
        printf("Measurement probes:\n");
        for(auto& p: probes) printf("< %g\t%g\t%g >\n", p.v0[0], p.v0[1], p.v0[2]);
    }
    
    
    ROOT::Math::Minimizer* myMin = NULL;        ///< minimizer core
    ROOT::Math::Functor f;
    vector<double> finalParams;
    
    /// set up fitter
    void init_fitter() {
        int nparam = myB.nparams() + 3*probes.size();
        printf("Initializing fitter for %i parameters\n", nparam);
        
        myMin = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
        //myMin = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "BFGS2");
        assert(myMin);
        f = ROOT::Math::Functor(this, &MapFitter::minFunc, nparam);
        
        myMin->SetMaxFunctionCalls(10000);  // for Minuit
        myMin->SetMaxIterations(1000);     // for GSL
        myMin->SetTolerance(0.01);
        myMin->SetFunction(f);
        
        set_init_params();
        
        const double* x = myMin->X();
        for(int i=0; i<nparam; i++) printf("\t%g",x[i]);
        printf("\n");
        myB.display(x);
    }
    
    void set_init_params() {
        for(int i=0; i<myB.nparams(); i++) {
            if(i==2) {
                if(offcenter) myMin->SetFixedVariable(i, "Bz", 368.);
                else myMin->SetVariable(i, "Bz", 368., 1.0);
            } else myMin->SetVariable(i, ("B"+to_str(i)).c_str(), 0, 0.1);
        }
        
        for(size_t n = 0; n < probes.size(); n++) {
            for(int i=0; i<3; i++) {
                if((!n && i==1) || (n==1 && i==0))  myMin->SetFixedVariable(3*n + i + myB.nparams(), ("p"+to_str(n)+"_"+to_str(i)).c_str(), probes[n].v0[i]);
                else myMin->SetVariable(3*n + i + myB.nparams(), ("p"+to_str(n)+"_"+to_str(i)).c_str(), probes[n].v0[i], 0.01);
            }
        }
    }
    
    /// fit field parameters
    void fit_field() {
        if(!myMin) init_fitter();
    
        printf("\nMinimizing...\n");
        myMin->Minimize();
        printf("\nDone!\n\n");
        
        const double* x = myMin->X();
        minFunc(x);
        finalParams.resize(myB.nparams() + 3*probes.size());
        for(int i=0; i<myB.nparams(); i++) finalParams[i] = x[i];
        myB.display(x);
    }
};

int main(int, char**) {
    
    string mapname = "1107";
    int nrot = 12; // 4 for 90 degree, 12 for 30 degree
    
    MapFitter MF;
    MF.offcenter = (nrot == 12);
    MF.set_probes();
    MF.setup_circle_rots(nrot);
    MF.init_fitter();
    
    string basedir = "/home/mpmendenhall/Documents/aCORN/FieldMaps/";
    string infname = basedir+"/MapData_Bmap"+mapname+".txt";
    FILE * datOut = fopen ((basedir+"Processed_"+mapname+".txt").c_str(),"w");
    assert(datOut);
    fprintf(datOut, "# Field map processed from '%s'\n", infname.c_str());
    fprintf(datOut, "# %zu probes at %zu angles per z\n", MF.probes.size(), MF.rots.size());
    fprintf(datOut, "# z[cm]\tt[s]\tT[C]\tBx\tBy\tBz");
    for(size_t i=0; i<MF.probes.size(); i++) fprintf(datOut, "\tP%zux\tP%zuy\tP%zuz", i, i, i);
    fprintf(datOut, "\n");
    
    auto datlines = split(loadFileString(infname),"\n");
    for(auto& l: datlines) {
        MF.set_probes();
        MF.set_init_params();
        MF.load_circle_data(l, nrot);
        MF.fit_field();
        MF.display_probes();
        printf("----------------------\n\n");
        
        fprintf(datOut,"%g\t%g\t%g", MF.scan_z, MF.scan_t, MF.scan_T);
        for(int i=0; i<3; i++) fprintf(datOut,"\t%g", MF.finalParams[i]);
        for(auto& p: MF.probes) for(int i=0; i<3; i++) fprintf(datOut,"\t%g", p.v0[i]);
        fprintf(datOut,"\n");
    }
    
    fclose(datOut);
    
    return EXIT_SUCCESS;
}

