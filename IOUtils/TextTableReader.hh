/// \file TextTableReader.hh Reader for tabbed text tables output by Igor
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include <vector>
using std::vector;
#include <string>
using std::string;

class TableData {
public:
    TableData(size_t nc = 0): ncols(nc) { }
    void addRows(const vector<double>& v);
    size_t ncols;
    size_t nrows() const { return dat.size()/ncols; }
    double operator()(size_t r, size_t c) const { return dat[r*ncols+c]; }
    vector<double> dat;
    void display() const { printf("Table %zux%zu\n", dat.size()/ncols, ncols); }
};

class TextTableReader {
public:
    TextTableReader() { }
    void loadFile(const string& fname);
    void loadText(const string& s);
    size_t mincols = 2;
    vector<TableData> tables;
    void display() const { for(auto const& t: tables) t.display(); }
};
