#!/usr/bin/python

from pyx import *
from pyx import style
from pyx import graph
from pyx.color import rgb, hsb
from pyx.graph.style import symbol
from LinFitter import *
import sqlite3
import os
try:
    from scipy import stats
except:
    stats = None
    
def get_result_types(curs):
    """Get list of analysis result types in DB"""
    curs.execute("SELECT rowid,name,descrip FROM named_object WHERE type = 'AnaResult'")
    return curs.fetchall()

def load_result_data(curs, num):
    """Load most recent data for specified analysis type"""
    sdat = {}
    for row in curs.execute("SELECT start_s,end_s,time,value,err FROM analysis_results WHERE type_id = %i ORDER BY time DESC"%num):
        sdat[(row[0], row[1])] = row[2:]
    return sdat
    
def load_seriesplot(curs, num, yscale):
    sdat = load_result_data(curs, num)
    gdat = [ [s[0]/10000, sdat[s][1]*yscale, sdat[s][2]*yscale] for s in sdat if s[0]/10000 == s[1]/10000 and sdat[s][1] is not None]
    gdat.sort()
    return gdat

class SeriesPlotSetup:
    def __init__(self, rt):
        self.num = rt[0]
        self.name = rt[1]
        self.scale = 1.
        self.ytitle = rt[2]
        self.basepath = os.environ["ACORN_SUMMARY"]+"/SeriesSummary/"
        self.gdat = []
        os.system("mkdir -p " + self.basepath)
        
        rescale = { "total_time":(1/3600.,"run time [h]"),
                "total_counts":(1e-6,"reduced events ($\\times 10^6$)")}
        if self.name in rescale:
            self.scale = rescale[self.name][0]
            self.ytitle = rescale[self.name][1]
         
    def make_gdat(self, curs):
        self.gdat = load_seriesplot(curs, self.num, self.scale)
        
    def make_graph(self, curs):
        
        self.xaxis = graph.axis.lin(title = "Series")
        self.yaxis = graph.axis.lin(title = self.ytitle)
        if not self.gdat:
            self.make_gdat(curs)
        g = graph.graphxy(width=15, height=6, x=self.xaxis, y=self.yaxis, key=graph.key.key(pos="tl"))
        g.plot(graph.data.points(self.gdat,x=1,y=2,dy=3,title=None),
               [graph.style.errorbar(), graph.style.symbol(symbol.circle,size=0.1,symbolattrs=[deco.filled([color.rgb(1,1,1)])])])
        
        self.auxfits(g)
        
        g.writetofile(self.basepath+self.name+".pdf")
        
    def summarize_dat(self):
        sx = 0
        sxw = 0
        sw = 0
        for x in self.gdat:
            sx += x[1]
            if not x[2]:
                x[2] = 1
            sxw += x[1]/x[2]**2
            sw += 1./x[2]**2
        print "Total =", sx, "Average =", sx/len(self.gdat), "Weighted =", sxw/sw, "N =", len(self.gdat)
        print
        
    def auxfits(self,g):
        if self.name in ["fid_asym","obs_asym"]:
            self.fit_line(g)
            
    def fit_line(self,g,nterms=1):
        LF = LinearFitter(terms=[polyterm(i) for i in range(nterms)])
        if self.gdat[0][2]:
            LF.fit(self.gdat, cols=(0,1,2), errorbarWeights = True)
        else:
            LF.fit(self.gdat, cols=(0,1))
        
        chi2 = LF.chisquared()
        nu = LF.nu()        
        ttl = "$"+LF.toLatex(cfmt=".3g")+"$; $\\chi^2/\\nu = %.1f/%i$"%(chi2,nu)
        if stats:
           ttl +=  " ($P = %.2f$)"%stats.chisqprob(chi2,nu)
           
        g.plot(graph.data.points(LF.fittedpoints(True), x=1,y=3,title=ttl),
               [graph.style.line(lineattrs=[rgb.red])])


def combo_series_dat(slist):
    csentries = {}
    for S in slist:
        for x in S.gdat:
            csentries.setdefault(x[0],[])
            csentries[x[0]] += x[1:]
    return [ [k,]+csentries[k] for k in csentries if len(csentries[k]) == 2*len(slist)]        

def calc_asym(r0,dr0,r1,dr1):
    s = r0 + r1
    ds = sqrt(dr0**2 + dr1**2)
    a = (r0-r1)/(r0+r1)
    da = 2*sqrt((r0*dr1)**2 + (r1*dr0)**2)/s**2
    return (s,ds,a,da)
    
def makeAsymmetrySeries(SP1,SP2):
    srate = SeriesPlotSetup((0,"fid_rate","wishbone fiducial rate [Hz]"))
    sasym = SeriesPlotSetup((0,"fid_asym","wishbone fiducial asymmetry"))
    gdat = combo_series_dat((SP1,SP2))
    for g in gdat:
       (s,ds,a,da) =  calc_asym(*g[1:])
       srate.gdat.append( [g[0],s,ds] )
       sasym.gdat.append( [g[0],a,da] )
    return srate, sasym

def plot_correlated(SP1, SP2):
    gdat = combo_series_dat((SP1,SP2))
    g = graph.graphxy(width=12,height=12, x=SP1.yaxis, y=SP2.yaxis, key=graph.key.key(pos="tl"))
    g.plot(graph.data.points(gdat,x=2,dx=3,y=4,dy=5,title=None),[graph.style.errorbar(), graph.style.symbol(symbol.circle,size=0.1,symbolattrs=[deco.filled([color.rgb(1,1,1)])])])
    
    LF = LinearFitter(terms=[polyterm(i) for i in range(2)])
    LF.fit(gdat, cols=(1,3,4), errorbarWeights = True)
    chi2 = LF.chisquared()
    nu = LF.nu()        
    ttl = "$"+LF.toLatex(cfmt=".3g")+"$; $\\chi^2/\\nu = %.1f/%i$"%(chi2,nu)
    if stats:
        ttl +=  " ($P = %.2f$)"%stats.chisqprob(chi2,nu)
    g.plot(graph.data.points(LF.fittedpoints(True), x=1,y=3,title=ttl), [graph.style.line(lineattrs=[rgb.red])])

    g.writetofile(SP1.basepath+SP2.name+"_VS_"+SP1.name+".pdf")
    
if __name__ == "__main__":
    conn = sqlite3.connect(os.environ["ACORN_DB"])
    curs = conn.cursor()
    
    rtplist = get_result_types(curs)
    rtpdat = {}
    for rt in rtplist:
        print rt
        rtpdat[rt[1]] = SeriesPlotSetup(rt)
        rtpdat[rt[1]].make_gdat(curs)
        
    rtpdat["fid_rate"],rtpdat["fid_asym"] = makeAsymmetrySeries(rtpdat["wb_fast_fiducial"], rtpdat["wb_slow_fiducial"])
    
    for k in rtpdat:
        print k
        rtpdat[k].make_graph(curs)
        rtpdat[k].summarize_dat()

    cplots = [("wb_bg","wb_fg"),
              ("wb_fast_fiducial", "wb_slow_fiducial"),
              ("fid_rate","obs_asym"),
              ("wb_bg","obs_asym"),
              ("fid_rate","wb_fg"),
              ("fid_asym","obs_asym")]
    for cplot in cplots:
        plot_correlated(rtpdat[cplot[0]], rtpdat[cplot[1]])
        