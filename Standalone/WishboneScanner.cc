#include "WishbonePlugin.hh"
#include "ReducedDataScanner.hh"
#include "PathUtils.hh"
#include "OutputManager.hh"
#include "strutils.hh"

int main(int, char**) {
    
    RunNum series = 999;
    
    ReducedDataScanner RDS;
    int nfailed = 0;
    int i = 0;
    while(nfailed < 5) {
        if(RDS.addRun(RunID(series,i))) nfailed = 0;
        else nfailed++;
        i++;
    }
    
    OutputManager OM("Wishbone", getEnvSafe("ACORN_WISHBONE")+"/Series_"+itos(series));
    WishboneAnalyzer WA(&OM);
    
    WA.loadProcessedData(RDS);
    WA.makeOutput();
    
    return 0;
}
