#!/usr/bin/python

from math import *
from LinFitter import *
from pyx import *
from pyx import style
from pyx import graph
from pyx.color import rgb, hsb
from pyx.graph.style import symbol
import os

def rainbow(n, b=1.0):
    return [ hsb((1.0*x)/n,1,b) for x in range(n) ]

class Gordon_Table:
    def __init__(self, fname, mul=1.0, basedir="/home/mpmendenhall/Documents/aCORN/Reference/CalculatedFieldmaps/"):
        self.mul = mul
        self.tbldat = [[float(x) for x in l.split(",")] for l in open(basedir+fname,"r").readlines()]
        self.nradii = len(self.tbldat[0])/3
        self.zs = [x[0] for x in self.tbldat]
        
    def get_radius_dat(self, r, col=2, mul = 1.0):
        return [ (x[0]*1e-2, x[3*r+col]*1e4*mul*self.mul) for x in self.tbldat]
    
    def add_field(self, GT, mul=1.0):
        assert len(self.tbldat) == len(GT.tbldat)
        for (n,l) in enumerate(self.tbldat):
            for m in range(len(l)-1):
                l[m+1] += mul*GT.tbldat[n][m+1]*GT.mul
        
class MPM_Table:
    def __init__(self, fname, basedir=os.environ["ROTSHIELD_OUT"]+"/aCORN/"):
        self.tbldat = [[float(x) for x in l.split()] for l in open(basedir+fname+"/Fieldmap.txt","r").readlines() if l[0] != "#"]
    def get_zdat(self):
        dz = 1.4981 # 1.4981 in nominal baseline
        return [(x[2]+dz,x[5]) for x in self.tbldat]
    
def fitdict(xydat):
    return (lambda x,m=dict(xydat): m[x])

def make_target(zlist):
    t = []
    B0 = 368
    for z in zlist:
        if 1.10 < z < 2.52:
            t.append((z, B0, 20. + (2.4 < z)*100))
        elif 2.70 < z < 3.00:
            l = (z-2.7)/(3-2.7)
            t.append((z, 250*(1-l)+90*l, 1e-6 + (2.8 < z)*0.3))  # 250, 70 for Ed
        else:
            t.append((z, 0, 1e-12))
    return t

def load_Gordon_tables(with_Ed = False):
    GT = [Gordon_Table("Bmap1_20.csv"),Gordon_Table("Bmap1_20.csv")]
    
    for i in range(21,26):
        fmapname = "Bmap%i.csv"%i
        if with_Ed and 22 <= i <= 24:
            fmapname = "Bmap%iEd.csv"%i
        mul = 1.0
        if with_Ed:
            if i == 24:
                mul = 9.5/29.2
            if i == 25:
                mul = 0
        GT.append(Gordon_Table(fmapname))
        GT[0].add_field(GT[-1], mul=mul)
  
    for i in range(21,26):
        mul = 5*(with_Ed and i==24)
        fmapname = "Bmap_AxialTrim%iEd.csv"%i
        GT.append(Gordon_Table(fmapname))
        GT[0].add_field(GT[-1], mul=mul)
        
    return GT

def transverse_from_axial(f, r, ds = 10):
    s0 = 0
    outdat = []
    while s0 + ds < len(f):
        dz = f[s0+ds][0] - f[s0][0]
        z0 = 0.5*(f[s0+ds][0] + f[s0][0])
        dB = f[s0+ds][1] - f[s0][1]
        s0 += ds
        outdat.append((z0, -0.5*r*dB/dz))
    return outdat

def transverse_field_study():
    GT = load_Gordon_tables()
    
    Brange = 0.2
    gMagF = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="position [m]",min=0,max=3),
        y=graph.axis.lin(title="$B_r$ [Gauss]", max=Brange, min=-Brange),
        key = graph.key.key(pos="tl", columns=2))
    
    nr = 7
    cls = rainbow(nr)
    for r in range(nr):
        ttl = "r = %i cm"%r
        gMagF.plot(graph.data.points(GT[0].get_radius_dat(r, 1),x=1,y=2,title=ttl),[graph.style.line([style.linewidth.Thick, cls[r]])])
        
    M = MPM_Table("Sym_Fat_2")
    Bt = transverse_from_axial(M.get_zdat(), 0.01*(nr-1), 1) #GT[0].get_radius_dat(0), 0.01*(nr-1))
    gMagF.plot(graph.data.points(Bt,x=1,y=2,title="MPM r=%i"%(nr-1)),[graph.style.line([style.linewidth.Thick])])
     
    gMagF.writetofile("MagF_transverse.pdf")

def field_shaper(with_Ed = False):
   
    GT = load_Gordon_tables(with_Ed)
    
    B0 = 370
    Bmin = 200
    if with_Ed:
        Bmin = 0
    
    gMagF = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="position [m]",min=0,max=3),
        y=graph.axis.lin(title="$B_z$ [Gauss]", min = Bmin, max = B0+10),
        key = graph.key.key(pos="bl"))
    gMagFNarrow = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="position [m]",min=0,max=3),
        y=graph.axis.lin(title="$B_z$ [Gauss]", min = 350, max = 380))
     
    if with_Ed:
        targf = make_target([p[0] for p in GT[0].get_radius_dat(0)])
        LF = LinearFitter(terms=[fitdict(G.get_radius_dat(0)) for G in GT[1:]], fixparams = {1:1.4, 3:1.4, 5:1.4, 8:5})
        LF.fit(targf, cols=(0,1,2))
        print LF.coeffs
        gMagF.plot(graph.data.points(LF.fittedpoints(),x=1,y=3,title=None),[graph.style.symbol(symbol.circle, size=0.02,)])
        gMagF.plot(graph.data.points(LF.fittedpoints(),x=1,y=2,title=None),[graph.style.line([style.linewidth.thin, style.linestyle.dotted])])
        
        gMagFNarrow.plot(graph.data.points(LF.fittedpoints(),x=1,y=3,title=None),[graph.style.symbol(symbol.circle, size=0.02,)])
        gMagFNarrow.plot(graph.data.function("y(x) = 368"),[graph.style.line([style.linewidth.thin, style.linestyle.dotted])])
         
    cls = rainbow(len(GT))
    for (n,G) in enumerate(GT):
        ttl = None
        if n and with_Ed:
            ttl = "Scaled by %g"%LF.coeffs[n-1]
        #elif with_Ed and not n:
        #    continue
        lstyle = graph.style.line([style.linewidth.Thick, cls[n]])
        if n > 6:
            lstyle = graph.style.line([style.linewidth.thick, style.linestyle.dotted, cls[n]])
        gMagF.plot(graph.data.points(G.get_radius_dat(0, mul = 1.0+(n>1)*0),x=1,y=2,title=ttl),[lstyle])
        if not n:
           gMagFNarrow.plot(graph.data.points(G.get_radius_dat(0, mul = 1.0+(n>1)*0),x=1,y=2,title=ttl),[graph.style.line([style.linewidth.Thick, cls[n]])])
        
    if not with_Ed:
        #M = MPM_Table("Baseline")
        M = MPM_Table("Baseline")
        M2 = MPM_Table("Sym_Fat_2")
        gMagF.plot(graph.data.points(M.get_zdat(),x=1,y=2,title=None),[graph.style.line([style.linewidth.thin])])
        gMagF.plot(graph.data.points(M2.get_zdat(),x=1,y=2,title=None),[graph.style.line([style.linewidth.thin, style.linestyle.dashed])])
    
    gMagF.writetofile("MagF.pdf")
    gMagFNarrow.writetofile("MagF_Detail.pdf")
    
if __name__=="__main__":
    #transverse_field_study()
    field_shaper(True)
