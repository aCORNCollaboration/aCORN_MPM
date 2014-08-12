#include "RunMetadata.hh"
#include "strutils.hh"

#include <fstream>
#include <iostream>
#include <cassert>

// metadata CSV version 6:
//
// 0    1           2               3    4                5              6              7     8         9
// Date,Julian Date,Wishbone Series,Tier,Start Run Number,End Run Number,Run Time (Sec),Slope,Intercept,Field
// 2/23/2013,13054,531,1,414,529,55661.9,0.0181865,41.364,-1

MetadataDB::MetadataDB(const std::string& fname) {
    if(!fname.size()) return;
    
    std::cout << "Loading metadata file '" << fname << "'...\n";
    
    std::ifstream infl(fname);
    
    bool pastHeader = false;
    while (infl) {
        std::string s;
        if (!getline(infl, s)) break;
        
        // check if we have gotten past the headers, ending at "Data,..." line
        if(!pastHeader) {
            pastHeader = (s[0] == 'D');
            continue;
        }
        
        std::vector<double> v = sToDoubles(s);
        if(v.size()<10) {
            std::cerr << "Bad line in metadata file: '" << s << "'\n";
            continue;
        }
        
        RunMetadata M;
        M.run.first = v[2];
        M.tier = RunMetadata::DataTier(v[3]);
        M.slope = v[7];
        M.intercept = v[8];
        M.field = RunMetadata::Field(v[9]);
        for(RunNum rn = v[4]; rn <= v[5]; rn++) {
            series[M.run.first].push_back(rn);
            M.run.second = rn;
            md[M.run] = M;
        }
    }
    
    if (!infl.eof())
        std::cerr << "Warning: EOF not reached in metadata file!\n";
    std::cout << "Loaded metadata for " << md.size() << " runs.\n";
}

RunMetadata MetadataDB::getRun(RunID r) const {
    std::map<RunID,RunMetadata>::const_iterator it = md.find(r);
    assert(it != md.end());
    if(it != md.end()) return it->second;
    return RunMetadata();
}
