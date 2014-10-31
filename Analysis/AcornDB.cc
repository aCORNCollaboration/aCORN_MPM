#include "AcornDB.hh"
#include "PathUtils.hh"
#include "SMExcept.hh"
#include "BaseDataScanner.hh"
#include <string.h>
#include <stdio.h>

void errorLogCallback(void*, int iErrCode, const char* zMsg){
    fprintf(stderr, "SQL error (%d): %s\n", iErrCode, zMsg);
}


AcornDB* AcornDB::myDB = NULL;

AcornDB& AcornDB::ADB() {
    if(!AcornDB::myDB) {
        sqlite3_config(SQLITE_CONFIG_LOG, &errorLogCallback, NULL);
        myDB = new AcornDB();
    }
    return *myDB;
}

AcornDB::AcornDB() {
    std::string dbname = getEnvSafe("ACORN_DB");
    printf("Opening SQLite3 DB '%s'...\n",dbname.c_str());
    int err = sqlite3_open(dbname.c_str(), &db);
    if(err) {
        SMExcept e("failed_db_open");
        e.insert("dbname",dbname);
        e.insert("message",sqlite3_errmsg(db));
        sqlite3_close(db);
        db = NULL;
        throw e;
    }
}

AcornDB::~AcornDB() {
    if(db) {
        for(auto it = statements.begin(); it != statements.end(); it++) sqlite3_finalize(*it);
        sqlite3_close(db);
    }
}

int AcornDB::setQuery(const char* qry, sqlite3_stmt*& stmt) {
    int rc = sqlite3_prepare_v2(db, qry, strlen(qry), &stmt, NULL);
    if(rc != SQLITE_OK) {
        SMExcept e("failed_query");
        e.insert("message",sqlite3_errmsg(db));
        throw(e);
    }
    return rc;
}

int combo_runid(RunID rn) { return 10000*rn.first + rn.second; }

void AcornDB::getPMTSumCal(RunID rn, double& slope, double& intercept) {
    
    static const char* qry = "SELECT slope,intercept FROM pmt_sum_cal WHERE start_s <= ?1 AND ?1 <= end_s ORDER BY end_s-start_s ASC";
    static sqlite3_stmt* stmt = NULL;
    if(!stmt) setQuery(qry, stmt);
    
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
    
    static const char* qry = "SELECT start_s,end_s FROM data_series WHERE ?2 >= start_s AND ?1 <= end_s AND type = 'wishbone' AND tier = ?3";
    static sqlite3_stmt* stmt = NULL;
    if(!stmt) setQuery(qry, stmt);
    
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

void AcornDB::loadPMTcal(RunID start, RunID end, int n, double sigPerPE, double sigPerMeV) {    
    static const char* qry = "INSERT OR REPLACE INTO pmt_gaincal(start_s, end_s, pmt, sigPerPE, sigPerMeV) VALUES (?1, ?2, ?3, ?4, ?5)";
    static sqlite3_stmt* stmt = NULL;
    if(!stmt) setQuery(qry, stmt);
    
    sqlite3_bind_int(stmt, 1, combo_runid(start));
    sqlite3_bind_int(stmt, 2, combo_runid(end));
    sqlite3_bind_int(stmt, 3, n);
    sqlite3_bind_double(stmt, 4, sigPerPE);
    sqlite3_bind_double(stmt, 5, sigPerMeV);
    
    sqlite3_step(stmt);
    sqlite3_reset(stmt);
}

void AcornDB::getPMTcal(RunID rn, vector<double>& sigPerPE, vector<double>& sigPerMeV) {
    static const char* qry = "SELECT pmt,sigPerPE,sigPerMeV FROM pmt_gaincal WHERE start_s <= ?1 AND ?1 <= end_s ORDER BY end_s-start_s ASC";
    static sqlite3_stmt* stmt = NULL;
    if(!stmt) setQuery(qry, stmt);

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
