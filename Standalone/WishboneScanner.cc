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

int main(int argc, char** argv) {
    aCORN_MPM::display_version();
    MPMUtils::display_version();
    
    if(argc != 2) {
        printf("Please supply a series number for analysis.\n");
        return EXIT_FAILURE;
    }

    setupSlideStyle();
    gROOT->ForceStyle();
    
    int series = atoi(argv[1]);
    
    // series 0: merge all wishbone data into run total
    if(series == 0) {
        printf("Merging series data...\n");
        OutputManager OM("Wishbone", getEnvSafe("ACORN_WISHBONE")+"/Combined/");
        string setname = argv[1];
        if(setname == "0") {
            WishboneAnalyzer WA(&OM,"Wishbone");
            WA.mergeDir(getEnvSafe("ACORN_WISHBONE"));
        } else {
            WishboneAnalyzer WA(&OM,"Wishbone_"+setname);
            vector<string> flist;
            for(auto s: AcornDB::ADB().groupSeries(setname))
                flist.push_back(getEnvSafe("ACORN_WISHBONE")+"/Series_"+to_str(s)+"/Series_"+to_str(s)+".root");
            WA.addFiles(flist);
            WA.makeOutput();
        }
        AcornDB::closeDB();
        return EXIT_SUCCESS;
    }
    
    string wbname = "Series_"+to_str(abs(series));
    OutputManager OM("NameUnused", getEnvSafe("ACORN_WISHBONE"));
    
    // negative series number: re-process output plots
    if(series < 0) {
        WishboneAnalyzer WA(&OM, wbname, OM.basePath+"/"+wbname+"/"+wbname+".root");
        WA.makeOutput();
        AcornDB::closeDB();
        return EXIT_SUCCESS;
    }
    
    ReducedDataScanner RDS(series >= 1519);
    if(!RDS.addRuns(AcornDB::ADB().seriesRuns(series))) {
        if(series > -3000) {
            printf("No wishbone runs in DB for series %u;",series);
            auto v = RDS.findSeriesRuns(series);
            if(!v.size()) {
                printf(" neither any runs on disk! Aaack!\n");
                return EXIT_FAILURE;
            }
            if(series == 3044) v.resize(420);
            RDS.addRuns(v);
            printf(" adding %zu runs found on disk.\n",v.size());
        } else {
            printf("Series %u contains no useful runs. Analysis stopped.\n", series);
            return EXIT_FAILURE;
        }
    }
    if(!RDS.nEvents) return EXIT_FAILURE;
    
    WishboneAnalyzer WA(&OM, wbname);

    WA.loadProcessedData(RDS);
    WA.makeOutput();
    
    AcornDB::closeDB();
    return EXIT_SUCCESS;
}
