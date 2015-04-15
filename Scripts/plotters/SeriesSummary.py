#!/usr/bin/python

from pyx import *
from pyx import style
from pyx import graph
from pyx.color import rgb, hsb
from pyx.graph.style import symbol
from LinFitter import *
import sqlite3
import os

def get_result_types(curs):
    """Get list of analysis result types in DB"""
    curs.execute("SELECT rowid,name,descrip FROM named_object WHERE type = 'AnaResult'")
    return curs.fetchall()

def load_result_data(curs, rt):
    """Load most recent data for specified analysis type"""
    sdat = {}
    for row in curs.execute("SELECT start_s,end_s,time,value,err FROM analysis_results WHERE type_id = %i ORDER BY time DESC"%rt[0]):
        sdat[(row[0], row[1])] = row[2:]
    return sdat
    
def plot_series_results(curs, rt, yinfo):
    """Generate series summary for analysis result"""
    basepath = os.environ["ACORN_SUMMARY"]+"/SeriesSummary/"
    os.system("mkdir -p " + basepath)
    
    sdat = load_result_data(curs, rt)
    gdat = [ [s[0]/10000, sdat[s][1]*yinfo[0], sdat[s][2]*yinfo[0]] for s in sdat if s[0]/10000 == s[1]/10000]
    
    g = graph.graphxy(width=15,height=6,
        x=graph.axis.lin(title = 'Series'),
        y=graph.axis.lin(title = yinfo[1]))
    
    g.plot(graph.data.points(gdat,x=1,y=2,dy=3),[graph.style.errorbar(), graph.style.symbol(symbol.circle,size=0.1,symbolattrs=[deco.filled([color.rgb(1,1,1)])])])
    g.writetofile(basepath+rt[1]+".pdf")
    
    sx = 0
    sxw = 0
    sw = 0
    for x in gdat:
        sx += x[1]
        if not x[2]:
            x[2] = 1
        sxw += x[1]/x[2]**2
        sw += 1./x[2]**2
    print "Total =", sx, "Average =", sx/len(gdat), "Weighted =", sxw/sw, "N =", len(gdat)
    print
    
if __name__ == "__main__":
    conn = sqlite3.connect(os.environ["ACORN_DB"])
    curs = conn.cursor()
    
    rescale = {"total_time":(1/3600.,"run time [h]")}
    
    for rt in get_result_types(curs):
        print rt
        plot_series_results(curs, rt, rescale.get(rt[1],(1.,rt[2])))