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
    def __init__(self, hdms, vdms, name=None, r=None):
        self.h = dms2r(*hdms)
        self.v = dms2r(*vdms)
        self.r = r
        self.name = name
        self.angular_cxn = (0,0)
    def offset_from(self,rpos):
        return ( (self.h-(rpos.h + rpos.angular_cxn[0]))*self.r, (self.v-(rpos.v + rpos.angular_cxn[1]))*self.r )

def theodolite_band_reading():
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
        
        
        
        
        
        
        
        
        
        
        
        

if __name__ == "__main__":
    
    # distances [inches] relative to bottom electron collimator
    rdists = {"bottom cage":    -14.,
              "bottom electron":0,
              "top electron":   14.750,
              "bottom mirror":  17.250,
              "top mirror":     34.625,
              "bottom proton":  40.000,
              "top proton":     95.000,
              "top cage":       120.000 }
    
    center_line = None
    
    if 0:
        r0 = 33.820
        center_line = ReticlePos((0,0,32),(90,0,28),1)
        rets = [ReticlePos((0,0,33),(90,0,23),"top electron"),
                ReticlePos((0,1,33),(90,0,28),"bottom mirror"),
                ReticlePos((0,0,54),(90,0,5),"top mirror"),
                ReticlePos((0,0,36),(90,0,32),"bottom proton"),
                ReticlePos((0,0,33),(90,0,31),"top proton")]
        
        center_line = ReticlePos((0,0,14),(90,0,32),1)
        rets = [ReticlePos((0,0,15),(90,0,41),"bottom electron"),
                ReticlePos((0,0,19),(90,0,35),"top electron"),
                ReticlePos((0,1,14),(90,0,36),"bottom mirror"),
                ReticlePos((0,0,33),(90,0,10),"top mirror"),
                ReticlePos((0,0,18),(90,0,41),"bottom proton"),
                ReticlePos((0,0,14),(90,0,36),"top proton")]
        
        # alignement in bore
        r0 = 54.7458
        center_line = ReticlePos((360,0,2),(89,56,24),1)
        rets = [ReticlePos((359,59,34),(89,56,32),"bottom electron"),
                ReticlePos((359,59,18),(89,57,0), "bottom mirror"),
                ReticlePos((359,59,38),(89,57,1), "top mirror"),
                ReticlePos((359,59,51),(89,56,42),"bottom proton"),
                ReticlePos((359,59,57),(89,56,36),"top proton")]
        
        # re-alignment in bore after removing bottom proton reticule
        r0 = 51.156
        rets = [ReticlePos((360,0,0),  (89,54,25),"bottom cage"),
                ReticlePos((359,58,59),(89,54,25),"bottom electron"),
                ReticlePos((359,59,0), (89,54,55),"bottom mirror"),
                ReticlePos((359,59,21),(89,54,57),"top mirror"),
                ReticlePos((359,59,39),(89,54,25),"top proton"),
                ReticlePos((359,59,59),(89,54,24),"top cage")]
        center_line = rets[-1]

        # re-alignment check MPM
        r0 = 51.156
        center_line = ReticlePos((0,43,39),(89,54,24),1)            # nominal centerline
        rets = [ReticlePos((0,43,37),(89,54,30),"bottom cage"),
                ReticlePos((0,42,37),(89,54,29),"bottom electron"),
                ReticlePos((0,42,33),(89,54,56),"bottom mirror"),
                ReticlePos((0,43,2), (89,54,56),"top mirror"),
                ReticlePos((0,43,16),(89,54,21),"top proton"),
                ReticlePos((0,43,38),(89,54,26),"top cage")]
        
        # re-alignment check MPM 20150131
        r0 = 50.704
        center_line = ReticlePos((0,4,15),(89,53,29),1)            # arb.
        rets = [ReticlePos((0,4,23),(89,53,15),"bottom cage"),
                ReticlePos((0,4,57),(89,53,54),"bottom electron"),
                ReticlePos((0,5,20),(89,54,4), "bottom mirror"),
                ReticlePos((0,5,2), (89,53,38),"top mirror"),
                ReticlePos((0,4,28),(89,53,47),"top proton"),
                ReticlePos((0,4,17),(89,53,27),"top cage")]
        
        # insert alignment 20150327 MPM
        #r0 = ???
        rets = [ReticlePos((0,0,4),    (89,55,8), "bottom cage"),
                ReticlePos((359,59,39),(89,55,32),"bottom electron"),
                ReticlePos((359,59,50),(89,55,59),"bottom mirror"),
                ReticlePos((0,0,7),    (89,55,44),"top mirror"),
                ReticlePos((359,59,57),(89,55,18),"top proton"),
                ReticlePos((0,0,0),    (89,55,5), "top cage")]
        center_line = rets[-1]
        
        # insert alignment 20150327 TH
        rets = [ReticlePos((0,0,3),    (89,55,6), "bottom cage"),
                ReticlePos((359,59,40),(89,55,29),"bottom electron"),
                ReticlePos((359,59,50),(89,55,57),"bottom mirror"),
                ReticlePos((0,0,4),    (89,55,43),"top mirror"),
                ReticlePos((359,59,53),(89,55,18),"top proton"),
                ReticlePos((0,0,0),    (89,55,7), "top cage")]
        center_line = rets[0]
        
        # insert alignment 20150608
        r0 = 48.102 - rdists["bottom cage"]
        rets = [ReticlePos((0,9,31), (89,53,47),"bottom cage"),
                ReticlePos((0,9,50), (89,54,43),"bottom electron"),
                ReticlePos((0,10,10),(89,54,42),"bottom mirror"),
                ReticlePos((0,10,0), (89,54,20),"top mirror"),
                ReticlePos((0,9,41), (89,54,12),"top proton"),
                ReticlePos((0,9,33), (89,53,57),"top cage")]
        center_line = rets[-1]
    
    if 0:
        # 20151216 in-bore alignment after opening up (TH); theodolite viewing from NNE
        r0 = 53.34 - rdists["bottom cage"]
        rets = [ReticlePos((0,0,0), (89,53,40),"bottom cage"),
                ReticlePos((0,0,0), (89,54,36),"bottom electron"),
                ReticlePos((0,0,28),(89,54,42),"bottom mirror"),
                ReticlePos((0,0,27),(89,54,23),"top mirror"),
                ReticlePos((0,0,5), (89,54,12),"top proton"),
                ReticlePos((0,0,3), (89,54,0), "top cage")]
        
    elif 0:
        # 20151216 re-check (MPM)
        r0 = 53.34 - rdists["bottom cage"]
        rets = [ReticlePos((0,0,12),(89,53,28),"bottom cage"),
                ReticlePos((0,0,13),(89,54,31),"bottom electron"),
                ReticlePos((0,0,38),(89,54,38),"bottom mirror"),
                ReticlePos((0,0,34),(89,54,15),"top mirror"),
                ReticlePos((0,0,21),(89,54,10),"top proton"),
                ReticlePos((0,0,14),(89,53,57),"top cage")]
        
    elif 1:
        # redone on more solid base (MPM)
        r0 = 53.34 - rdists["bottom cage"]
        rets = [ReticlePos((0,14,45),(89,53,15),"bottom cage"),
                ReticlePos((0,15,13),(89,54,22),"bottom electron"),
                ReticlePos((0,15,58),(89,54,34),"bottom mirror"),
                ReticlePos((0,16,14),(89,54,8), "top mirror"),
                ReticlePos((0,16,25),(89,54,4), "top proton"),
                ReticlePos((0,16,24),(89,53,51),"top cage")]
        
    #elif 0:
        # 20160112 aligning insert outside bore
    
    # assign distances
    for r in rets:
        if r.name in rdists:
            r.r = rdists[r.name] + r0
    
    ############################################
    # re-set centerline along bearing cage axis!
    if not center_line:
        center_line = rets[-1]
    def realign_center(cl, r0, r1):
        """re-align centerline through reticules"""
        cl.angular_cxn = [ (r1.offset_from(cl)[i]-r0.offset_from(cl)[i])/(r1.r-r0.r) for i in range(2) ]
    realign_center(center_line, rets[0], rets[-1])
    
    g = graph.graphxy(width=10,height=10,
        x=graph.axis.lin(title="horizontal offset [mils]", min=-20, max=20),
        y=graph.axis.lin(title="vertical offset [mils]", min=-5, max=35),
        key = graph.key.key(pos="tc",columns=2))
    
    # symbols for drawing points
    rsymbs = {"bottom cage": symbol.diamond,
              "bottom electron": symbol.triangle,
              "top electron": symbol.diamond, # not used with bottom cage
              "bottom mirror": symbol.plus,
              "top mirror": symbol.cross,
              "bottom proton": symbol.square, # not used with top cage
              "top proton": symbol.circle,
              "top cage": symbol.square }
    
    ret0 = rets[0].offset_from(center_line)
    for (n,r) in enumerate(rets):
        gdat = r.offset_from(center_line)
        gdat = [ (1000*(gdat[0]-ret0[0]), 1000*(gdat[1]-ret0[1])) ] # inches to mils
        symb = rsymbs.get(r.name, None)
        ttl = r.name + ": (%.1f, %.1f)"%gdat[0]
        g.plot(graph.data.points(gdat, x=1, y=2, title=ttl),[graph.style.symbol(symb)])
    
    if 0:
        # compare previous set    
        r0 = 51.156
        center_line = ReticlePos((0,43,39),(89,54,24),1)            # nominal centerline
        rets = [ReticlePos((0,43,37),(89,54,30),"bottom cage"),
                ReticlePos((0,42,37),(89,54,29),"bottom electron"),
                ReticlePos((0,42,33),(89,54,56),"bottom mirror"),
                ReticlePos((0,43,2), (89,54,56),"top mirror"),
                ReticlePos((0,43,16),(89,54,21),"top proton"),
                ReticlePos((0,43,38),(89,54,26),"top cage")]
        names = ["bottom cage", "bottom electron", "bottom mirror", "top mirror",  "top proton",   "top cage"]
        symbs = [symbol.diamond, symbol.triangle,  symbol.plus,    symbol.cross, symbol.circle, symbol.square]
        for (n,r) in enumerate(rets):
            gdat = [rotate(r.offset_from(center_line),-3*pi/4)]
            symb = symbs[n]
            g.plot(graph.data.points(gdat, x=1, y=2, title=None),[graph.style.symbol(symb,symbolattrs=[rgb.red])])
        
        
    g.writetofile("/home/mpmendenhall/Desktop/crosshairs_20150608.pdf")
    