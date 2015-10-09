/// \file TextTableReader.cc

#include "TextTableReader.hh"
#include "StringManip.hh"
#include <cassert>

void TableData::addRows(const vector<double>& v) { 
    assert(!(v.size()%ncols));
    dat.insert(dat.end(), v.begin(), v.end());
}

void TextTableReader::loadText(const string& s) {
    auto lines = split(s,"\r\n");
    for(auto const& l: lines) {
        auto vals = sToDoubles(l);
        if(vals.size() < mincols) continue;
        if(!tables.size() || tables.back().ncols != vals.size()) tables.push_back(TableData(vals.size()));
        tables.back().addRows(vals);
    }
}

void TextTableReader::loadFile(const string& fname) {
    loadText(loadFileString(fname));
}
