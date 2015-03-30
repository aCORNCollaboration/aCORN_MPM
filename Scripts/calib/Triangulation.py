#!/usr/bin/python

from math import *

from pyx import *
from pyx import style
from pyx import graph
from pyx.color import rgb, hsb
from pyx.graph.style import symbol

def d2dms(d):
    """degrees to degrees, minutes, seconds"""
    if d<0:
        d += 360
    return (floor(d), floor(60*(d%1)), 60*((60*d)%1))

def dms2d(d,m,s):
    """degrees, minutes, seconds to degrees"""
    a = d + m/60. + s/3600.
    if a > 180:
        a = a-360.
    return a

def d2r(d):
    """degrees to radians"""
    return d*pi/180.

def dms2r(d,m,s):
    """degrees, minutes, seconds to radians"""
    return d2r(dms2d(d,m,s))
    
class PerpPos:
    """Single-theodolite location relative to perpendicular plane"""
    def __init__(self, a1, a2, dy):
        self.x0 = dy * cos(a2)*sin(pi/2+a1)/sin(a2-a1)
        self.yref = self.ypos(a1)
        print "Theodolite at x =",self.x0," above band top at y =",self.yref
    def ypos(self, a):
        """y position of point on plane at angle a"""
        return -tan(a)*self.x0

class Triangulator:
    """Two-theodolite, single-plane triangulator"""
    def __init__(self, PP0, PP1):
        self.x0 = PP0.x0
        self.x1 = PP1.x0
        self.dy = PP0.yref - PP1.yref
        print "Theodolite pair spaced dy =", self.dy
    def xypos(self,a1,a2):
        """(x,y) position of point for theodolite angles a1,a2"""
        # -y/(x0+x) = tan(a1)
        # (dy-y)/(x1+x) = tan(a2)
        t1,t2 = tan(a1),tan(a2)
        x = (self.dy+t1*self.x0-t2*self.x1)/(t2-t1)
        y = -t1*(self.x0+x)
        return (x,y)
    def angles(self,x,y):
        """(x,y) position to theodolite angles"""
        return (atan2(-y,self.x0+x), atan2(self.dy-y,self.x1+x))

def offset_dist(p1,p2,r):
    return sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)*r

def rotate(xy,theta):
    return (xy[0]*cos(theta) - xy[1]*sin(theta), xy[0]*sin(theta) + xy[1]*cos(theta))
    
    
class ReticlePos:
    def __init__(self,hdms,vdms,r):
        self.h = dms2r(*hdms)
        self.v = dms2r(*vdms)
        self.r = r
    def offset_from(self,rpos):
        return ( (self.h-rpos.h)*self.r, (self.v-rpos.v)*self.r )
        
if __name__ == "__main__":
    if 0:
        w_band = 1./2.54    # width of copper band, inches
        
        # theodolite readings at top and bottom of band
        t0top = dms2r(93,4,5) - pi/2
        t0bot = dms2r(93,54,5) - pi/2
        t1top = dms2r(135,4,5) - pi/2
        t1bot = dms2r(135,54,5) - pi/2
        
        PP0 = PerpPos(t0top, t0bot, w_band)
        PP1 = PerpPos(t1top, t1bot, w_band)
        T = Triangulator(PP0, PP1)
        
        t0p = dms2r(94,5,6) - pi/2
        t1p = dms2r(135,12,8) - pi/2
        (x,y) = T.xypos(t0p,t1p)
        (a1,a2) = T.angles(x,y)
        print "angles",(t0p,t1p),"-> position",(x,y),"-> angles",(a1,a2)

    r_BC = -14000.
    r_BEC = 0.
    r_TEC = 14750.
    r_BM = 17250.
    r_TM = 34625.
    r_BPC = 40000.
    r_TPC = 95000.
    r_TC = 120000.
    
    r0 = 33820
    center_line = ReticlePos((0,0,32),(90,0,28),1)
    rets = [ReticlePos((0,0,33),(90,0,23),r0+14750),
            ReticlePos((0,1,33),(90,0,28),r0+17250),
            ReticlePos((0,0,54),(90,0,5),r0+34625),
            ReticlePos((0,0,36),(90,0,32),r0+40000),
            ReticlePos((0,0,33),(90,0,31),r0+95000)]
    
    center_line = ReticlePos((0,0,14),(90,0,32),1)
    rets = [ReticlePos((0,0,15),(90,0,41),r0),          # bottom EC
            ReticlePos((0,0,19),(90,0,35),r0+14750),    # top EC
            ReticlePos((0,1,14),(90,0,36),r0+17250),    # bottom mirror
            ReticlePos((0,0,33),(90,0,10),r0+34625),    # top mirror
            ReticlePos((0,0,18),(90,0,41),r0+40000),    # bottom proton
            ReticlePos((0,0,14),(90,0,36),r0+95000)]    # top proton
    names = ["bottom electron", "top electron","bottom mirror","top mirror","bottom proton","top proton"]
    symbs = [symbol.triangle, symbol.circle, symbol.plus, symbol.cross, symbol.square, symbol.diamond]
    
    # alignement in bore
    r0 = 54745.8
    center_line = ReticlePos((360,0,2),(89,56,24),1)            # top bearing cage
    rets = [ReticlePos((359,59,34),(89,56,32),r0),              # bottom E
            ReticlePos((359,59,18),(89,57,0),r0+17250),         # bottom EM
            ReticlePos((359,59,38),(89,57,1),r0+34625),         # top EM
            ReticlePos((359,59,51),(89,56,42),r0+40000),        # bottom P
            ReticlePos((359,59,57),(89,56,36),r0+95000)]        # top P
    names = ["bottom electron", "bottom mirror", "top mirror", "bottom proton", "top proton"]
    symbs = [symbol.triangle, symbol.plus, symbol.cross, symbol.square, symbol.diamond]
    
    # re-alignment in bore after removing bottom proton reticule
    r0 = 51156.
    center_line = ReticlePos((360,0,0),(89,54,24),1)            # ???
    rets = [ReticlePos((360,0,0),(89,54,25),r0+r_BC),
            ReticlePos((359,58,59),(89,54,25),r0+r_BEC),
            ReticlePos((359,59,0),(89,54,55),r0+r_BM),
            ReticlePos((359,59,21),(89,54,57),r0+r_TM),
            ReticlePos((359,59,39),(89,54,25),r0+r_TPC),
            ReticlePos((359,59,59),(89,54,24),r0+r_TC)]
    names = ["bottom cage", "bottom electron", "bottom mirror", "top mirror",  "top proton",   "top cage"]
    symbs = [symbol.diamond, symbol.triangle,  symbol.plus,    symbol.cross, symbol.circle, symbol.square]

    # re-alignment check MPM
    r0 = 51156.
    center_line = ReticlePos((0,43,39),(89,54,24),1)            # nominal centerline
    rets = [ReticlePos((0,43,37),(89,54,30),r0+r_BC),
            ReticlePos((0,42,37),(89,54,29),r0+r_BEC),
            ReticlePos((0,42,33),(89,54,56),r0+r_BM),
            ReticlePos((0,43,2), (89,54,56),r0+r_TM),
            ReticlePos((0,43,16),(89,54,21),r0+r_TPC),
            ReticlePos((0,43,38),(89,54,26),r0+r_TC)]
    
    # re-alignment check MPM 20150131
    r0 = 50704.
    center_line = ReticlePos((0,4,15),(89,53,29),1)            # arb.
    rets = [ReticlePos((0,4,23),(89,53,15),r0+r_BC),
            ReticlePos((0,4,57),(89,53,54),r0+r_BEC),
            ReticlePos((0,5,20),(89,54,4),r0+r_BM),
            ReticlePos((0,5,2), (89,53,38),r0+r_TM),
            ReticlePos((0,4,28),(89,53,47),r0+r_TPC),
            ReticlePos((0,4,17),(89,53,27),r0+r_TC)]
    
    # insert alignment 20150327 MPM
    #r0 = ???
    center_line = ReticlePos((0,0,0),(89,55,5),1)
    rets = [ReticlePos((0,0,4),(89,55,8),r0+r_BC),
            ReticlePos((359,59,39),(89,55,32),r0+r_BEC),
            ReticlePos((359,59,50),(89,55,59),r0+r_BM),
            ReticlePos((0,0,7), (89,55,44),r0+r_TM),
            ReticlePos((359,59,57),(89,55,18),r0+r_TPC),
            ReticlePos((0,0,0),(89,55,05),r0+r_TC)]
    
    # insert alignment 20150327 TH
    center_line = ReticlePos((0,0,0),(89,55,6),1)
    rets = [ReticlePos((0,0,3),(89,55,6),r0+r_BC),
            ReticlePos((359,59,40),(89,55,29),r0+r_BEC),
            ReticlePos((359,59,50),(89,55,57),r0+r_BM),
            ReticlePos((0,0,4), (89,55,43),r0+r_TM),
            ReticlePos((359,59,53),(89,55,18),r0+r_TPC),
            ReticlePos((0,0,0),(89,55,07),r0+r_TC)]
    
    g = graph.graphxy(width=10,height=10,
        x=graph.axis.lin(title="horizontal offset [mils]", min=-20, max=20),
        y=graph.axis.lin(title="vertical offset [mils]", min=-10, max=30),
        key = graph.key.key(pos="tc",columns=2))
    
    for (n,r) in enumerate(rets):
        gdat = [r.offset_from(center_line)]
        symb = symbs[n]
        ttl = names[n] + ": (%.1f, %.1f)"%gdat[0]
        g.plot(graph.data.points(gdat, x=1, y=2, title=ttl),[graph.style.symbol(symb)])
    
    if 0:
        # compare previous set    
        r0 = 51156.
        center_line = ReticlePos((0,43,39),(89,54,24),1)            # nominal centerline
        rets = [ReticlePos((0,43,37),(89,54,30),r0+r_BC),
                ReticlePos((0,42,37),(89,54,29),r0+r_BEC),
                ReticlePos((0,42,33),(89,54,56),r0+r_BM),
                ReticlePos((0,43,2), (89,54,56),r0+r_TM),
                ReticlePos((0,43,16),(89,54,21),r0+r_TPC),
                ReticlePos((0,43,38),(89,54,26),r0+r_TC)]
        names = ["bottom cage", "bottom electron", "bottom mirror", "top mirror",  "top proton",   "top cage"]
        symbs = [symbol.diamond, symbol.triangle,  symbol.plus,    symbol.cross, symbol.circle, symbol.square]
        for (n,r) in enumerate(rets):
            gdat = [rotate(r.offset_from(center_line),-3*pi/4)]
            symb = symbs[n]
            g.plot(graph.data.points(gdat, x=1, y=2, title=None),[graph.style.symbol(symb,symbolattrs=[rgb.red])])
        
        
    g.writetofile("/home/mpmendenhall/Desktop/crosshairs.pdf")
    