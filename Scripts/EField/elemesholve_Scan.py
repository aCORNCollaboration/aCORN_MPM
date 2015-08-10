#!/usr/bin/python

from math import *
from pyx import *
from pyx import style
from pyx import graph
from pyx.color import rgb, hsb
from pyx.graph.style import symbol

def rainbow(n, b=1.0):
    return [ hsb((1.0*x)/n,1,b) for x in range(n) ]
    
def rainbowDict(keys, b=1.0):
    n = len(keys)
    knew = [k for k in keys]
    knew.sort()
    return dict([ (k,hsb((1.0*x)/n,1,b)) for (x,k) in enumerate(knew) ])
        
class GridPt:
    def __init__(self, fline):
        d = [float(x) for x in fline.split()]
        self.x = d[:3]
        self.phi = d[3]
        self.E = d[4:]
    def dist(self, gp):
        return sqrt(sum([(self.x[i]-gp.x[i])**2 for i in range(3)]))
        
def load_scans(fname):
    
    scans = [ ]
    for p in [GridPt(l) for l in open(fname).readlines()]:
        if scans and p.dist(scans[-1][-1]) < 0.2:
            scans[-1].append(p)
        else:
            scans.append([p,])
            
    return scans


    
def plot_scans(fname):
    scans = load_scans(fname)
    cols = rainbow(len(scans))
    
    gPhi = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="position [cm]"),
        y=graph.axis.log(title="potential [kV]"))
    
    gEt = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="position [cm]"),
        y=graph.axis.lin(title="transverse field [kV/cm]", min = -0.01, max = 0.01))
    
    gEa = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="position [cm]", min = -1, max = 1),
        y=graph.axis.lin(title="axial field [kV/cm]"))
    
    for (n,sc) in enumerate(scans):
        gdat = [(p.x[2], p.phi, sqrt(p.E[0]**2 + p.E[1]**2), p.E[0], p.E[1], p.E[2]) for p in sc]
        
        plotl = [graph.style.line([cols[n]])]
        plotl2 = [graph.style.line([cols[n], style.linestyle.dotted])]
        
        gPhi.plot(graph.data.points(gdat,x=1,y=2), plotl)
        gEt.plot(graph.data.points(gdat,x=1,y=4), plotl)
        #gEt.plot(graph.data.points(gdat,x=1,y=5), plotl2)
        gEa.plot(graph.data.points(gdat,x=1,y=6), plotl)
        
    gPhi.writetofile("FieldScans.pdf")
    gEt.writetofile("FieldScans_Et.pdf")
    gEa.writetofile("FieldScans_Ea.pdf")
    
if __name__ == "__main__":
    plot_scans("~/Applications/elemesholve-bld/scan.txt")
    