#!/usr/bin/python

from sys import path
path.append("../plotters")

from CalDB import *
from LinFitter import *
from pyx import *
from pyx import style
from pyx import graph
from pyx.color import rgb, hsb
from pyx.graph.style import symbol


class SourcePeak:
    def __init__(self):
        pass
    
def collect_sdata(s0, s1):
    conn = sqlite3.connect(os.environ["ACORN_DB"])
    curs = conn.cursor()
    #                    0       1       2       3       4      5       6        7       8
    curs.execute("SELECT series, pktype, dttype, center, sigma, height, dcenter, dsigma, dheight FROM calib_peaks WHERE %i <= series AND series <= %i"%(s0,s1))
    pkdat = curs.fetchall()
    conn.close()
    
    sdata = {}
    for p in pkdat:
        pk = SourcePeak()
        pk.center = p[3]
        pk.sigma = p[4]
        pk.height = p[5]
        pk.dcenter = p[6]
        pk.dsigma = p[7]
        pk.dheight = p[8]
        sdata.setdefault((p[0],p[1]),{})[p[2]] = pk
        
    return sdata
    
def make_sourcecal(sdata):
    g = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="PMT sum [ADC channels]", min=0),
        y=graph.axis.lin(title="Expected energy [keV]", min=0),
        key = graph.key.key(pos="tl"))
     
    gdat = [(pks["ADC_sum"].center, pks["ADC_sum"].dcenter, pks["MC"].center) for pks in sdata.itervalues()]
    g.plot(graph.data.points(gdat, x=1, dx=2, y=3, title=None),[graph.style.symbol(symbol.circle)])
    
    LF = LinearFitter(terms=[polyterm(i) for i in range(2)])
    LF.fit(gdat,cols=(0,2))
    g.plot(graph.data.points(LF.fitcurve(0,4e4), x=1, y=2, title="$%s$"%LF.toLatex(cfmt=".5g")),[graph.style.line(lineattrs=[style.linewidth.Thick])])
    
    g.writetofile("ADCSum_Cal.pdf")
    
if __name__ == "__main__":
    make_sourcecal(collect_sdata(1327,1329));
    
#UPDATE pmt_sum_cal SET slope = 0.02529, intercept = 29.41 WHERE start_s = 13250000 AND end_s = 13309999;