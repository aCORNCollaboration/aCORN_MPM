#include "WishbonePlugin.hh"
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
    
    // series 0: merge all wishbone data into run total
    if(series == 0) {
        OutputManager OM("Wishbone", dropLast("/"+strip(getEnvSafe("ACORN_WISHBONE"),"/"),"/"));
        WishboneAnalyzer WA(&OM,"Wishbone");
        printf("Merging series data in '%s'...\n", WA.basePath.c_str());
        WA.mergeDir();
        return 0;
    }
    
    OutputManager OM("Wishbone", getEnvSafe("ACORN_WISHBONE"));
    
    // negative series number: re-process output plots
    if(series < 0) {
        series = -series;
        string wbname = "Series_"+to_str(series);
        WishboneAnalyzer WA(&OM, wbname, OM.basePath+"/"+wbname+"/"+wbname);
        WA.makeOutput();
        return 0;
    }
    
    ReducedDataScanner RDS(series >= 1519);
    if(!RDS.addRuns(AcornDB::ADB().seriesRuns(series))) {
        printf("Series %u contains no useful runs. Analysis stopped.\n", series);
        return 0;
    }
    
    WishboneAnalyzer WA(&OM, "/Series_"+to_str(series));
        
    WA.loadProcessedData(RDS);
    WA.makeOutput();
    
    return 0;
}
