from math import *
from scipy import special
import numpy

from pyx import *
from pyx import style
from pyx import graph
from pyx.color import rgb, hsb
from pyx.graph.style import symbol


def I0sum(x, nmax = 60):
    s = 0
    for n in range(nmax):
        s += x**(2*n)/(4**n * gamma(n+1)**2)
    return s

def f0(n):
    return -(-1)**n/(4*pi*n)


class MirrorField:
    def __init__(self):
        self.V0 = 1
        self.L = 0.7
        self.w = 6.8/7.0
        self.r0 = 5.5
    def f1(self,n):
        return     (cos(pi*n*self.w)-(-1)**n)/(2*pi*n)
    def f2(self,n):
        return (self.w*cos(pi*n*self.w)-(-1)**n)/(4*pi*n) - sin(pi*n*self.w)/(2*pi*n)**2
    def kn(self,n):
        return 2*pi*n/self.L
    def I0scaled(self,r,n):
        return special.iv(0, self.kn(n)*r)
    def cn_sharp(self,n):
        return -4 * self.V0 * (f0(n)-0.5*self.f1(n)) / self.I0scaled(self.r0,n)
        #return -4 * self.V0 * (f0(n)-0.5*self.f1(n))
    def cn_smooth(self,n):
        ff2 = self.f2(n)
        return -4*self.V0 * (f0(n)-ff2+self.w/(1-self.w)*(0.5*self.f1(n)-ff2)) / self.I0scaled(self.r0,n)
        #return -4*self.V0 * (f0(n)-ff2+self.w/(1-self.w)*(0.5*self.f1(n)-ff2))

    def eval_sharp(self, r, z, nt = 10):
        s = 0
        for n in range(nt)[1:]:
            s += self.cn_sharp(n) * sin(self.kn(n)*z) * self.I0scaled(r,n)
        return s
    
    def eval_smooth(self, r, z, nt = 10):
        s = 0
        for n in range(nt)[1:]:
            s += self.cn_smooth(n) * sin(self.kn(n)*z) * self.I0scaled(r,n)
        return s

if __name__ == "__main__":
    MF = MirrorField()
    for n in range(20)[1:]:
        #print "\t%i\t& %f\t& %f\t& %f\t& %g\t& %g\t& %g \\\\"%(n, f0(n), MF.f1(n), MF.f2(n), MF.I0scaled(MF.r0,n), MF.cn_sharp(n), MF.cn_smooth(n))
        print "\t%i\t& %f\t& %f\t& %g\t\\\\"%(n, MF.cn_sharp(n), MF.cn_smooth(n), MF.I0scaled(MF.r0,n))

    g = graph.graphxy(width=10,height=6,
            x=graph.axis.lin(title="z position [cm]"),
            y=graph.axis.lin(title="potential $V^*(z,r)/V_0$"))
    
    for r in [5.5, 5.4, 5.3, 5.2, 5.1, 5.0]:
        gdat = [ (z, MF.eval_smooth(r, z), MF.eval_sharp(r, z)) for z in numpy.linspace(-1.4,1.4,400) ]
        #g.plot(graph.data.points(gdat,x=1,y=2),[graph.style.line(lineattrs=[style.linewidth.Thick,rgb.red])])
        g.plot(graph.data.points(gdat,x=1,y=3),[graph.style.line(lineattrs=[style.linewidth.thin,rgb.blue])])
        g.writetofile("mfield.pdf")

    g2 = graph.graphxy(width=10,height=6,
            x=graph.axis.lin(title="distance from edge $r_0-r$ [cm]"),
            y=graph.axis.log(title="potential $|V^*(0.25, r)|/V_0$")) #, max = 0.5, min = 1e-6))    
    gdat = [ (r, -MF.eval_sharp(5.5-r, 0.25)) for r in numpy.linspace(0,2,200) ]
    g2.plot(graph.data.points(gdat,x=1,y=2),[graph.style.line(lineattrs=[style.linewidth.thin,rgb.blue])])
    g2.writetofile("mradial.pdf")    

