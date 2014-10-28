#!/usr/bin/python

import sqlite3
import os

def load_metadata_csv(curs, csv_filename):
    """Import csv "metadata" file into DB"""
    
    curs.execute("DELETE FROM wishbone_series WHERE 1")
    curs.execute("DELETE FROM pmt_sum_cal WHERE 1")
    
    # 0    1           2               3    4                5              6              7     8         9
    # Date,Julian Date,Wishbone Series,Tier,Start Run Number,End Run Number,Run Time (Sec),Slope,Intercept,Field
    past_header = False
    for l in [l.strip().split(",") for l in open(csv_filename,"r").readlines()]:
        # scan to actual data listing
        if not past_header:
            past_header = l[0]=="Date"
            continue
        
        s = int(l[2])
        r0 = int(l[4])
        r1 = int(l[5])
        s0 = 10000*s+r0
        s1 = 10000*s+r1
        t = l[0].split('/')
        t = "%s-%02i-%02i"%(t[2],int(t[0]),int(t[1]))
        
        cmd = "INSERT INTO wishbone_series VALUES (%i, %i, strftime('%%s','%s'), %s, %s, %s)"%(s0,s1,t,l[6],l[3],l[9])
        print cmd
        curs.execute(cmd)
        curs.execute("INSERT INTO pmt_sum_cal VALUES(%i,%i,%s,%s)"%(s0,s1,l[7],l[8]))
    
    
    
if __name__ == "__main__":
    conn = sqlite3.connect(os.environ["ACORN_DB"])
    curs = conn.cursor()
    
    load_metadata_csv(curs,os.environ["ACORN_METADATA"]);
    
    conn.commit()
    conn.close()
