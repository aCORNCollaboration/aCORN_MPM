from math import *
from scipy import special
import numpy

from pyx import *
from pyx import style
from pyx import graph
from pyx.color import rgb, hsb
from pyx.graph.style import symbol

# pre-compute bessel function zeros
bessel_j0n = [0,]+list(special.jn_zeros(0, 1000))

class BesselCalcs:
    """Field for flat mismatched-voltage endcap on cylinder"""
    def __init__(self, r0 = 1):
        # radius scale units
        self.r0 = r0
        # pre-compute coefficients
        self.ncoeffs = 101
        self.cns = [0 for n in range(self.ncoeffs)]
        
    def add_endcircle(self, r, V):
        """Add circular patch of given radius and voltage"""
        for n in range(self.ncoeffs)[1:]:
            self.cns[n] += V*self.cn(n, r)
        
    def cn(self,n,ri):
        """Expansion coefficient"""
        if self.r0 == ri:
            return 2/(bessel_j0n[n]*special.jv(1,bessel_j0n[n]))
        r = ri/self.r0
        return 2*special.jv(1,r*bessel_j0n[n])/(bessel_j0n[n]*special.jv(1,bessel_j0n[n])**2)
    
    def V(self, r, z):
        """Potential as a function of r,z position"""
        r /= self.r0
        z /= self.r0
        s = 0
        for (n,c) in enumerate(self.cns):
            if n==0:
                continue
            s += self.cns[n]*exp(-bessel_j0n[n]*z)*special.jv(0,bessel_j0n[n]*r)
        return s

    def Er(self, r, z):
        """Radial field component dV/dr"""
        r /= self.r0
        z /= self.r0
        s = 0
        for (n,c) in enumerate(self.cns):
            if n==0:
                continue
            s += self.cns[n]*exp(-bessel_j0n[n]*z)*bessel_j0n[n]*special.jvp(0,bessel_j0n[n]*r)
        return s/self.r0
    
    def Ez(self, r, z):
        """axial field component dV/dz"""
        r /= self.r0
        z /= self.r0
        s = 0
        for (n,c) in enumerate(self.cns):
            if n==0:
                continue
            s += -bessel_j0n[n]*self.cns[n]*exp(-bessel_j0n[n]*z)*special.jv(0,bessel_j0n[n]*r)
        return s/self.r0 

class TwoEndedBessel:
    def __init__(self,dz,BC1,BC2=None):
        if not BC2:
            BC2 = BesselCalcs(BC1.r0)   # grounded end
        for n in range(len(BC1.cns))[1:]:
            ez = exp(-bessel_j0n[n]*dz/BC1.r0)
            k = 1./(1-ez**2)
            BC1.cns[n], BC2.cns[n] = (BC1.cns[n]-ez*BC2.cns[n])*k, (-ez*BC1.cns[n]+BC2.cns[n])*k
        self.BC1 = BC1
        self.BC2 = BC2
        self.dz = dz
        
    def V(self, r, z):
        return self.BC1.V(r,z) + self.BC2.V(r, self.dz-z)
    def Er(self, r, z):
        return self.BC1.Er(r,z) + self.BC2.Er(r, self.dz-z)
    def Ez(self, r, z):
        return self.BC1.Ez(r,z) - self.BC2.Ez(r, self.dz-z)
    
if __name__ == "__main__":
    BC = BesselCalcs()
    print BC.j0n
    for n in range(200)[1:]:
        print "\t%i\t& %g\t\\\\"%(n, BC.cn(n))

    g = graph.graphxy(width=10,height=6,
            x=graph.axis.lin(title="radial position $r/r_0$"),
            y=graph.axis.lin(title="endcap contribution $V(r,z)/V_1$"))
    
    for z in [0.1*n for n in range(25)]:
        npts = 50
        if z==0:
            npts = 600
        gdat = [ (r, BC.V(r,z)) for r in numpy.linspace(0,1,400) ]
        g.plot(graph.data.points(gdat,x=1,y=2),[graph.style.line(lineattrs=[style.linewidth.thin,rgb.blue])])
        g.writetofile("endfield.pdf")

