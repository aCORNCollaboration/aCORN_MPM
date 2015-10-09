#!/usr/bin/python

import os
from math import *
from EndField import *

from pyx import *
from pyx import style
from pyx import graph
from pyx.color import rgb, hsb
from pyx.graph.style import symbol
import cPickle
from numpy import *

def rainbow(n, b=1.0):
    return [ hsb((1.0*x)/n,1,b) for x in range(n) ]
    
def rainbowDict(keys, b=1.0):
    n = len(keys)
    knew = [k for k in keys]
    knew.sort()
    return dict([ (k,hsb((1.0*x)/n,1,b)) for (x,k) in enumerate(knew) ])
        
class GridPt:
    """Grid scan point in field map"""
    def __init__(self, fline = "", mode="mpm"):
        if fline.strip():
            if mode == "mpm":
                self.from_mpmline(fline)
            elif mode == "brian":
                self.from_brianline(fline)
                
    def from_mpmline(self, fline):
        d = [float(x) for x in fline.split()]
        self.x = d[:3]
        self.phi = d[3]
        self.E = d[4:]
        self.init_calculated()
            
    def from_brianline(self, fline):
        d = [float(x) for x in fline.split()]
        self.x = d[:3]
        self.phi = 0
        self.E = d[3:]
        self.init_calculated()
        
    def init_calculated(self):
        self.Er = 0
        r = self.r()
        if r:
            self.Er = (self.E[0]*self.x[0] + self.E[1]*self.x[1])/r
        
    def dist(self, gp):
        """Distance from another point"""
        return sqrt(sum([(self.x[i]-gp.x[i])**2 for i in range(3)]))
    def r(self):
        """Radial coordinate"""
        return sqrt(self.x[0]**2 + self.x[1]**2)
    def add(self, pt):
        """add another point's field"""
        self.phi += pt.phi
        self.Er += pt.Er
        for i in range(3):
            self.E[i] += pt.E[i]
    def mul(self, c):
        """multiply field by a constant"""
        self.phi *= c
        self.Er *= c
        for i in range(3):
            self.E[i] *= c
    
class ScanLine:
    def __init__(self, pts):
        self.pts = pts
        self.npts = len(self.pts)
        # coordinate mean and standard deviation
        self.mu = [sum([p.x[i] for p in self.pts])/self.npts for i in range(3)]
        self.rms = [sqrt(sum([(p.x[i]-self.mu[i])**2 for p in self.pts])/self.npts) for i in range(3)]
        # identify varying axis
        for i in range(3):
            if self.rms[i] > 0.01:
                self.varaxis = i
    def r(self):
        return sqrt(self.mu[0]**2 + self.mu[1]**2)
            
def average_lines(slines):
    lens = [s.npts for s in slines]
    lens.sort()
    medlen = lens[len(lens)/2]
    print len(slines), medlen
    
    nsum = 0
    ss = None
    for s in slines:
        if s.npts != medlen:
            continue
        nsum += 1
        if not ss:
            ss = s
            continue
        for (n,p) in enumerate(s.pts):
            ss.pts[n].add(p)
    for p in ss.pts:
        p.mul(1./nsum)
    
    return ss
    
def load_scans(fname):
    """Load scan data, separated into lines"""
    
    # load separate scan lines
    scans = [ ]
    for p in [GridPt(l) for l in open(fname).readlines()]:
        if scans and p.dist(scans[-1][-1]) < 0.2:
            scans[-1].append(p)
        else:
            scans.append([p,])
    print "Found", len(scans), "scans."
    
    # group by radius
    scangroups = [ ]
    for s in scans:
        l = ScanLine(s)
        if scangroups and abs(l.r() - scangroups[-1][-1].r()) < 0.01:
            scangroups[-1].append(l)
        else:
            scangroups.append([l,])
    
    print "\tin", len(scangroups), "radius groups."
    return [average_lines(g) for g in scangroups]


    
def plot_scans(fname):
    scans = load_scans(fname)
    cols = rainbow(len(scans))
    
    gPhi = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="position [cm]", min = -0.5, max = 3),
        y=graph.axis.lin(title="potential [kV]", max = 0.04))
    
    gEt = graph.graphxy(width=12,height=8,
        x=graph.axis.lin(title="position [cm]", min = -8, max = 8),
        y=graph.axis.lin(title="transverse field [V/cm]", min = -4, max = 0.5),
        key=graph.key.key(pos="bl"))
    
    gEa = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="position [cm]", min = -1, max = 1),
        y=graph.axis.lin(title="axial field [V/cm]"))
    
    V = 4.1
    BC1 = BesselCalcs(5.5)
    BC1.add_endcircle(4.14, V)
    
    BC2 = BesselCalcs(10.)
    BC2.add_endcircle(10, 3*V)
    BC2.add_endcircle(6.5, -3*V)
    BC2.add_endcircle(3.928, V)
    BC2 = TwoEndedBessel(10, BC2)
    
    def bccombo(r,z):
        if z < 0:
            return BC1.Er(r,-z)
        return BC2.Er(r,z)
    
    for (n,sc) in enumerate(scans):
        print sc.r(), sc.npts
        gdat = [(p.x[2], p.phi, p.Er, p.E[0], p.E[1], p.E[2]) for p in sc.pts]
        gtitle = "$r = %.1f$ cm"%sc.r()
        
        plotl = [graph.style.line([cols[n], style.linewidth.thin])]
        plotl2 = [graph.style.line([cols[n], style.linestyle.dashed, style.linewidth.THin])]
        
        gPhi.plot(graph.data.points(gdat,x=1,y=2), plotl)
        gEt.plot(graph.data.points(gdat,x=1,y=3,title=gtitle), plotl)
        #gEt.plot(graph.data.points(gdat,x=1,y=5), plotl2)
        gEa.plot(graph.data.points(gdat,x=1,y=6), plotl)
        
        gdat_endsim = [(p.x[2], bccombo(p.r(),p.x[2])) for p in sc.pts[::5]]
        gEt.plot(graph.data.points(gdat_endsim,x=1,y=2,title=None), plotl2)
        
    gPhi.writetofile("FieldScans.pdf")
    gEt.writetofile("FieldScans_Et.pdf")
    gEa.writetofile("FieldScans_Ea.pdf")
    
if __name__ == "__main__":
    plot_scans(os.environ["APP_DIR"]+"/elemesholve-bld/scan.txt")

    
