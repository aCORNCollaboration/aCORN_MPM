/// \file AcornDB.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#include "AcornDB.hh"
#include "PathUtils.hh"
#include "SMExcept.hh"
#include "BaseDataScanner.hh"
#include <string.h>
#include <stdio.h>
#include <time.h>

AcornDB* AcornDB::myDB = NULL;

AcornDB& AcornDB::ADB() {
    if(!myDB) myDB = new AcornDB();
    return *myDB;
}

void AcornDB::closeDB() {
    if(myDB) {
        delete myDB;
        myDB = NULL;
    }
}

AcornDB::AcornDB(): SQLite_Helper(getEnvSafe("ACORN_DB")) {
}

int combo_runid(RunID rn) { return 10000*rn.first + rn.second; }

void AcornDB::getPMTSumCal(RunID rn, double& slope, double& intercept) {
    static sqlite3_stmt* stmt = loadStatement("SELECT slope,intercept FROM pmt_sum_cal WHERE start_s <= ?1 AND ?1 <= end_s ORDER BY end_s-start_s ASC");
    
    sqlite3_bind_int(stmt, 1, combo_runid(rn));
    int rc1 = sqlite3_step(stmt);
    if(rc1 == SQLITE_ROW) {
        slope = sqlite3_column_double(stmt, 0);
        intercept = sqlite3_column_double(stmt, 1);
    } else {
        SMExcept e("missing_pmtsumcal");
        e.insert("message",sqlite3_errmsg(db));
        throw(e);
    }
    sqlite3_reset(stmt);
}

vector<RunID> AcornDB::seriesRuns(RunNum S, DataTier T) {
    static sqlite3_stmt* stmt = loadStatement("SELECT start_s,end_s FROM data_series WHERE ?2 >= start_s AND ?1 <= end_s AND type = 'wishbone' AND tier = ?3");
    
    vector<RunID> v;
    
    sqlite3_bind_int(stmt, 1, 10000*S);
    sqlite3_bind_int(stmt, 2, 10000*(S+1)-1);
    sqlite3_bind_int(stmt, 3, T);
    while(sqlite3_step(stmt) == SQLITE_ROW) {
        int s0 = sqlite3_column_int(stmt, 0);
        int s1 = sqlite3_column_int(stmt, 1);
        for(int i=s0; i<=s1; i++) {
            RunID rn(i/10000, i%10000);
            if(rn.first == S) v.push_back(rn);
        }
    }
    sqlite3_reset(stmt);
    
    return v;
}

vector<RunNum> AcornDB::groupSeries(const string& sname) {
    static sqlite3_stmt* stmt = loadStatement("SELECT start_s,end_s FROM wishbone_groups WHERE name = ?1");
    sqlite3_bind_text(stmt, 1, sname.c_str(), -1, SQLITE_STATIC);
    int rc1 = sqlite3_step(stmt);
    vector<RunNum> v;
    if(rc1 == SQLITE_ROW) {
        RunNum s0 = sqlite3_column_int(stmt, 0);
        RunNum s1 = sqlite3_column_int(stmt, 1);
        static sqlite3_stmt* stmt2 = loadStatement("SELECT DISTINCT start_s/10000 FROM data_series WHERE start_s/10000 = end_s/10000 AND ?1 <= start_s/10000 AND start_s/10000 <= ?2 AND tier = 1;");
        sqlite3_bind_int(stmt2, 1, s0);
        sqlite3_bind_int(stmt2, 2, s1);
        while(sqlite3_step(stmt2) == SQLITE_ROW) v.push_back(sqlite3_column_int(stmt2, 0));
        sqlite3_reset(stmt2);
    }
    sqlite3_reset(stmt);
    return v;
}

void AcornDB::loadPMTcal(RunID start, RunID end, int n, double sigPerPE, double sigPerMeV) {    
    static sqlite3_stmt* stmt = loadStatement("INSERT OR REPLACE INTO pmt_gaincal(start_s, end_s, pmt, sigPerPE, sigPerMeV) VALUES (?1, ?2, ?3, ?4, ?5)");
    
    sqlite3_bind_int(stmt, 1, combo_runid(start));
    sqlite3_bind_int(stmt, 2, combo_runid(end));
    sqlite3_bind_int(stmt, 3, n);
    sqlite3_bind_double(stmt, 4, sigPerPE);
    sqlite3_bind_double(stmt, 5, sigPerMeV);
    
    sqlite3_step(stmt);
    sqlite3_reset(stmt);
}

void AcornDB::getPMTcal(RunID rn, vector<double>& sigPerPE, vector<double>& sigPerMeV) {
    static sqlite3_stmt* stmt = loadStatement("SELECT pmt,sigPerPE,sigPerMeV FROM pmt_gaincal WHERE start_s <= ?1 AND ?1 <= end_s ORDER BY end_s-start_s ASC");

    sigPerPE = vector<double>(N_E_PMT,0);
    sigPerMeV = vector<double>(N_E_PMT,0);
    sqlite3_bind_int(stmt, 1, combo_runid(rn));
    while(sqlite3_step(stmt) == SQLITE_ROW) {
        size_t n = sqlite3_column_int(stmt, 0);
        if(n >= N_E_PMT) continue;
        if(sigPerPE[n]) break;
        sigPerPE[n] = sqlite3_column_double(stmt, 1);
        sigPerMeV[n] = sqlite3_column_double(stmt, 2);
    }
    sqlite3_reset(stmt);
}

TGraphErrors* AcornDB::getRecal(RunID rn) {
    static sqlite3_stmt* stmt = loadStatement("SELECT graph_ID FROM wishbone_recal WHERE start_s <= ?1 AND ?1 <= end_s ORDER BY end_s-start_s ASC");
    
    sqlite3_bind_int(stmt, 1, combo_runid(rn));
    sqlite3_int64 gid = 0;
    if(sqlite3_step(stmt) == SQLITE_ROW) gid = sqlite3_column_int64(stmt, 0);
    sqlite3_reset(stmt);
    
    if(!gid) return NULL;
    return getGraph(gid);
}

void AcornDB::uploadPeak(const CalPeak& pk, bool replace) {
    static sqlite3_stmt* dstmt = loadStatement("DELETE FROM calib_peaks WHERE series = ?1 AND pktype = ?2 AND dttype = ?3");
    static sqlite3_stmt* istmt = loadStatement("INSERT INTO calib_peaks(series, pktype, dttype, center, sigma, height, dcenter, dsigma, dheight) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)");
    
    if(replace) {
        sqlite3_bind_int(dstmt, 1, pk.series);
        sqlite3_bind_text(dstmt, 2, pk.pktype.c_str(), -1, SQLITE_STATIC);
        sqlite3_bind_text(dstmt, 3, pk.dttype.c_str(), -1, SQLITE_STATIC);    
        sqlite3_step(dstmt);
        sqlite3_reset(dstmt);
    }
    
    sqlite3_bind_int(istmt, 1, pk.series);
    sqlite3_bind_text(istmt, 2, pk.pktype.c_str(), -1, SQLITE_STATIC);
    sqlite3_bind_text(istmt, 3, pk.dttype.c_str(), -1, SQLITE_STATIC); 
    sqlite3_bind_double(istmt, 4, pk.center.x);
    sqlite3_bind_double(istmt, 5, pk.sigma.x);
    sqlite3_bind_double(istmt, 6, pk.height.x);
    sqlite3_bind_double(istmt, 7, pk.center.err);
    sqlite3_bind_double(istmt, 8, pk.sigma.err);
    sqlite3_bind_double(istmt, 9, pk.height.err);
    sqlite3_step(istmt);
    sqlite3_reset(istmt);
}

TGraphErrors* AcornDB::getGraph(sqlite3_int64 gID) {
    static sqlite3_stmt* stmt = loadStatement("SELECT x,dx,y,dy FROM graph_points WHERE graph_id = ?1");
    
    TGraphErrors* g = new TGraphErrors();
    
    sqlite3_bind_int64(stmt, 1, gID);
    int idx = 0;
    while(sqlite3_step(stmt) == SQLITE_ROW) {
        g->SetPoint(idx, sqlite3_column_double(stmt, 0), sqlite3_column_double(stmt, 2));
        g->SetPointError(idx, sqlite3_column_double(stmt, 1), sqlite3_column_double(stmt, 3));
        idx++;
    }
    sqlite3_reset(stmt);

    return g;
}

sqlite3_int64 AcornDB::createNamed(const string& tp, const string& name, const string& descrip) {
    static sqlite3_stmt* stmt = loadStatement("INSERT INTO named_object(type, name, descrip) VALUES (?1, ?2, ?3)");
    
    sqlite3_bind_text(stmt, 1, tp.c_str(), -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 2, name.c_str(), -1, SQLITE_STATIC); 
    sqlite3_bind_text(stmt, 3, descrip.c_str(), -1, SQLITE_STATIC); 
    sqlite3_step(stmt);
    sqlite3_reset(stmt);
    
    return sqlite3_last_insert_rowid(db);
}

sqlite3_int64 AcornDB::uploadGraph(const TGraph* g, const string& name, const string& descrip) {
    if(!g) return 0;
    sqlite3_int64 gid = createNamed("graph",name,descrip);
    printf("Uploading graph '%s' [%s] with ID %lli\n", name.c_str(), descrip.c_str(), gid);
    
    static sqlite3_stmt* stmt = loadStatement("INSERT INTO graph_points(graph_id,x,dx,y,dy) VALUES (?1, ?2, ?3, ?4, ?5)");

    double x,y,dx,dy;
    dx = dy = 0;
    const TGraphErrors* ge = dynamic_cast<const TGraphErrors*>(g);
    beginTransaction();
    for(int i=0; i<g->GetN(); i++) {
        g->GetPoint(i,x,y);
        if(ge) { 
            dx = ge->GetErrorX(i);
            dy = ge->GetErrorY(i);
        }
        sqlite3_bind_int64(stmt, 1, gid);
        sqlite3_bind_double(stmt, 2, x);
        sqlite3_bind_double(stmt, 3, dx);
        sqlite3_bind_double(stmt, 4, y);
        sqlite3_bind_double(stmt, 5, dy);
        sqlite3_step(stmt);
        sqlite3_reset(stmt);
    }
    endTransaction();
    
    return gid;
}

sqlite3_int64 AcornDB::uploadRecal(const TGraph* g, RunID r0, RunID r1) {
    static sqlite3_stmt* stmt = loadStatement("INSERT INTO wishbone_recal(start_s, end_s, graph_ID) VALUES (?1, ?2, ?3)");
    
    sqlite3_int64 gid = uploadGraph(g, "energyRecal", "Energy re-calibration curve for " + to_str(r0) + " -- " + to_str(r1));
    sqlite3_bind_int(stmt, 1, combo_runid(r0));
    sqlite3_bind_int(stmt, 2, combo_runid(r1));
    sqlite3_bind_int64(stmt, 3, gid); 
    sqlite3_step(stmt);
    sqlite3_reset(stmt);
    
    return gid;
}

sqlite3_int64 AcornDB::getAnaResultType(const string& name, const string& descrip) {
    static sqlite3_stmt* stmt = loadStatement("SELECT rowid FROM named_object WHERE type = 'AnaResult' AND name = ?1");
    
    sqlite3_bind_text(stmt, 1, name.c_str(), -1, SQLITE_STATIC);
    sqlite3_int64 id = 0;
    if(sqlite3_step(stmt) == SQLITE_ROW) id = sqlite3_column_int64(stmt, 0);
    else id = createNamed("AnaResult", name, descrip);
    sqlite3_reset(stmt);
    return id;
}

void AcornDB::uploadAnaResult(sqlite3_int64 type_id, AnaResult R) {
    static sqlite3_stmt* stmt = loadStatement("INSERT INTO analysis_results(type_id, start_s, end_s, time, value, err) VALUES (?1, ?2, ?3, ?4, ?5, ?6)");
    
    sqlite3_bind_int64(stmt, 1, type_id);
    sqlite3_bind_int(stmt, 2, combo_runid(R.start));
    sqlite3_bind_int(stmt, 3, combo_runid(R.end));
    sqlite3_bind_int64(stmt, 4, R.time? R.time : time(NULL));
    sqlite3_bind_double(stmt, 5, R.value);
    sqlite3_bind_double(stmt, 6, R.err);
    sqlite3_step(stmt);
    sqlite3_reset(stmt);
}

void AcornDB::uploadAnaResult(const string& name, const string& descrip, AnaResult R) {
    printf("%s [%s]:\t%g ~ %g\n", descrip.c_str(), name.c_str(), R.value, R.err);
    uploadAnaResult(getAnaResultType(name,descrip), R);
}
