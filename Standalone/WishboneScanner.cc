/// \file WishboneScanner.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "WishbonePlugin.hh"
#include "ReducedDataScanner.hh"
#include "PathUtils.hh"
#include "OutputManager.hh"
#include "AcornDB.hh"
#include "StringManip.hh"
#include "GraphicsUtils.hh"
#include "CompileVersion.hh"
#include "aCornCompileVersion.hh"

#include <stdio.h>
#include <TStyle.h>
#include <TROOT.h>

class WishboneAnalyzer: public RunAccumulator {
public:
    WishboneAnalyzer(OutputManager* pnt, const std::string& nm = "Wishbone", const std::string& inflname = ""):
    RunAccumulator(pnt, nm, inflname) {
        myBuilders["WishbonePlugin"] = &myWishbonePluginBuilder;
        buildPlugins();
    }
    
    WishbonePluginBuilder myWishbonePluginBuilder;
};

int main(int argc, char** argv) {
    aCORN_MPM::display_version();
    MPMUtils::display_version();
    
    if(argc != 2) {
        printf("Please supply a series number for analysis.\n");
        return 0;
    }
    int series = atoi(argv[1]);
    
    setupSlideStyle();
    gROOT->ForceStyle();
    
    // series 0: merge all wishbone data into run total
    if(series == 0) {
        OutputManager OM("Wishbone", dropLast("/"+strip(getEnvSafe("ACORN_WISHBONE"),"/"),"/"));
        WishboneAnalyzer WA(&OM,"Wishbone");
        printf("Merging series data in '%s'...\n", WA.basePath.c_str());
        WA.mergeDir();
        return 0;
    }
    
    string wbname = "Series_"+to_str(abs(series));
    OutputManager OM("NameUnused", getEnvSafe("ACORN_WISHBONE"));
    
    // negative series number: re-process output plots
    if(series < 0) {
        WishboneAnalyzer WA(&OM, wbname, OM.basePath+"/"+wbname+"/"+wbname+".root");
        WA.makeOutput();
        return 0;
    }
    
    ReducedDataScanner RDS(series >= 1519);
    if(!RDS.addRuns(AcornDB::ADB().seriesRuns(series))) {
        if(series > -3000) {
            auto v = RDS.findSeriesRuns(series);
            if(series == 3044) v.resize(420);
            RDS.addRuns(v);
            printf("No wishbone runs in DB for series %u; adding %zu runs found on disk.\n",series,v.size());
        } else {
            printf("Series %u contains no useful runs. Analysis stopped.\n", series);
            return 0;
        }
    }
    
    WishboneAnalyzer WA(&OM, wbname);

    WA.loadProcessedData(RDS);
    WA.makeOutput();
    
    return 0;
}
