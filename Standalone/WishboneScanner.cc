#include "WishbonePlugin.hh"
#include "ReducedDataScanner.hh"
#include "PathUtils.hh"
#include "OutputManager.hh"
#include "RunMetadata.hh"
#include "strutils.hh"

#include <stdio.h>

int main(int argc, char** argv) {
    
    if(argc != 2) {
        printf("Please supply a series number for analysis.\n");
        return 0;
    }
    
    RunNum series = atoi(argv[1]);
    
    if(series == 0) {
        printf("Merging series data...\n");
        OutputManager OM("Wishbone", "/home/mpmendenhall/Data/");
        WishboneAnalyzer WA(&OM,"aCORN_Wishbone");
        WA.mergeDir();
        return 0;
    }
    
    ReducedDataScanner RDS(series >= 1519);
    if(!RDS.addRuns(MetadataDB::MDB.seriesRuns(series))) {
        printf("Series %u contains no useful runs. Analysis stopped.\n", series);
        return 0;
    }
    
    OutputManager OM("Wishbone", getEnvSafe("ACORN_WISHBONE"));
    WishboneAnalyzer WA(&OM, "/Series_"+itos(series));
    
    WA.loadProcessedData(RDS);
    WA.makeOutput();
    
    return 0;
}
