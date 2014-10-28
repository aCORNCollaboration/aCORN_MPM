#!/usr/bin/python

import sqlite3
import os

def load_metadata_csv(curs, csv_filename):
    """Import csv "metadata" file into DB"""
    
    curs.execute("DELETE FROM data_series WHERE type = 'wishbone'")
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
        
        cmd = "INSERT INTO data_series(start_s,end_s,type,date,runtime,tier,field) VALUES (%i, %i, 'wishbone', strftime('%%s','%s'), %s, %s, %s)"%(s0,s1,t,l[6],l[3],l[9])
        print cmd
        curs.execute(cmd)
        curs.execute("INSERT INTO pmt_sum_cal VALUES(%i,%i,%s,%s)"%(s0,s1,l[7],l[8]))
        
def load_cal_csv(curs, csv_filename):
    """Import calibration data spreadsheet into DB"""
    curs.execute("DELETE FROM data_series WHERE type != 'wishbone'")
    
    prevdate = ''
    
    #0 1    2      3        4       5      6         7     8      9         10          11     12     13    14    15               16    17
    # ,Date,Series,Strt Run,End Run,Source,Run Time ,Slope,dSlope,Intercept,dIntercept ,ThresE,ThresV,VoltE,VoltV,Lower Beta Trim ,Field,Comments
    for l in [[l.strip() for l in l.strip().split(",")] for l in open(csv_filename,"r").readlines()]:
        if len(l) < 17:
                continue
        if len(l[1].split('.'))==3:
                prevdate = l[1]
        else:
                l[1] = prevdate
        if not (l[2] and l[3] and l[4] and l[5]):
            print "Skip line",l
            continue
        if l[5][:4] == "Wish" or l[5] == "Source":
            continue
        l[5] = l[5].replace("(High Stat)","").strip()
        if l[5] not in ["Bi","Sn","Ce","Ba","Co","Cd","Bg","Bg*","Beam OFF WB","Reactor OFF WB"]:
            print "Unknown type '"+l[5]+"'"
                
        s = int(l[2].replace("s",""))
        r0 = int(l[3])
        try:
            r1 = int(l[4])
        except:
            r1 = 0
        s0 = 10000*s+r0
        s1 = 10000*s+r1
        
        insert_fields = {"start_s":"%i"%s0, "end_s":"%i"%s1,"type":"'%s'"%l[5]}
        
        t = l[1].split('.')
        t = "%s-%02i-%02i"%(t[2],int(t[0]),int(t[1]))
        insert_fields["date"] = "strftime('%%s','%s')"%t
        
        if l[11]:
            insert_fields["energy_thresh"] = "'%s'"%l[11]
        if l[12]:
            insert_fields["veto_thresh"] = "%g"%float(l[12])
        if l[13]:
            insert_fields["energy_volt"] = "'%s'"%l[13]
        if l[14]:
            insert_fields["veto_volt"] = "'%s'"%l[14]
        if l[15]:
            insert_fields["beta_trim"] = "%g"%float(l[15])
        if l[16]=="DN":
            insert_fields["field"] = "-1"
        if l[16]=="UP":
            insert_fields["field"] = "1"
        if len(l) >= 18 and l[17]:
            insert_fields["comments"] = "'%s'"%l[17]

        cmd = "INSERT INTO data_series("
        vkeys = insert_fields.keys()
        vkeys.sort()
        for k in vkeys:
            cmd += "%s,"%k
        cmd = cmd[:-1]+") VALUES ("
        for k in vkeys:
            cmd += "%s,"%insert_fields[k]
        cmd = cmd[:-1]+")"
        
        print cmd
        curs.execute(cmd)
    
if __name__ == "__main__":
    conn = sqlite3.connect(os.environ["ACORN_DB"])
    curs = conn.cursor()
    
    #load_metadata_csv(curs,os.environ["ACORN_METADATA"]);
    load_cal_csv(curs,"/home/mpmendenhall/Documents/aCORN/Database_for_Daily_Calibration_Full DataBase.csv");
    
    conn.commit()
    conn.close()
