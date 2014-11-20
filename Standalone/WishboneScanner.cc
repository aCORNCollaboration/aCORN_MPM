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
    
    if(series == 0) {
        printf("Merging series data...\n");
        OutputManager OM("Wishbone", "/home/mpmendenhall/Data/");
        WishboneAnalyzer WA(&OM,"aCORN_Wishbone");
        WA.mergeDir();
        return 0;
    }
    
    OutputManager OM("Wishbone", getEnvSafe("ACORN_WISHBONE"));
    
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
