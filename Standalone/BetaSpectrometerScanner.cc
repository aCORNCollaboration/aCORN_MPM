#include "PMTsPlugin.hh"
#include "ReducedDataScanner.hh"
#include "PathUtils.hh"
#include "OutputManager.hh"
#include "AcornDB.hh"
#include "StringManip.hh"

#include <stdio.h>
#include <TStyle.h>

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
