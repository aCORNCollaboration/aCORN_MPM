/// \file BetaSpectrometerScanner.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "PMTsPlugin.hh"
#include "ReducedDataScanner.hh"
#include "PathUtils.hh"
#include "OutputManager.hh"
#include "AcornDB.hh"
#include "StringManip.hh"

#include <stdio.h>
#include <TStyle.h>

class PMTsAnalyzer: public RunAccumulator {
public:
    PMTsAnalyzer(OutputManager* pnt, const std::string& nm = "PMTs", const std::string& inflname = ""):
    RunAccumulator(pnt, nm, inflname), myPMTsPluginBuilder(this) {
        myBuilders["PMTsPlugin"] = &myPMTsPluginBuilder;
        buildPlugins();
    }
    
    PMTsPluginBuilder myPMTsPluginBuilder;
};

int main(int argc, char** argv) {
    
    if(argc != 2) {
        printf("Please supply a series number for analysis.\n");
        return 0;
    }
    
    int series = atoi(argv[1]);
    
    gStyle->SetOptStat("");
    
    OutputManager OM("BetaSpectrometer", getEnvSafe("ACORN_WISHBONE"));
    
    ReducedDataScanner RDS(series >= 1519);
    auto v = RDS.findSeriesRuns(series);
    if(!RDS.addRuns(v)) return 0;
    
    PMTsAnalyzer PA(&OM, "/Series_"+to_str(series));
    
    PA.loadProcessedData(RDS);
    PA.makeOutput();
    
    return 0;
}
