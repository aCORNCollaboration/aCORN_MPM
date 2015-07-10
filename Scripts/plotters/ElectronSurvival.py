#!/usr/bin/python

from math import *
from LinFitter import *
from pyx import *
from pyx import style
from pyx import graph
from pyx.color import rgb, hsb
from pyx.graph.style import symbol

# l: collimator diameter
# d: distance of origin point from collimator center
# a: Larmor radius

def nintegrate(f, x1, x2, npts=20):
    """Simpson's Rule crude numerical integrator"""
    s = 0;
    for i in range(npts+1):
        x = x1 + (x2-x1)*float(i)/float(npts);
        w = 2.;
        if i%2:
            w = 4.
        if i==0 or i==npts:
            w = 1.
        s += f(x)*w;
    return s*(x2-x1)/(3.*npts);
        
def plane_survival_fraction(l, d, a):
    """Proportion of trajectories in plane that pass collimator"""
    if d >= l or 2*a >= d+l:
        return 0.
    if d + 2*a <= l:
        return 1.
    if d == 0:
        return 2*a < l
    x = (l**2 - d**2 - 2*l*a)/(2*d*a)
    if x <= -1:
        return 0
    if x >= 1:
       return 1
    return 1 - acos(x)/pi;

def sphere_survival_fraction(l, d, a0):
    """Proportion of trajectories emitted into sphere with max Larmor radius a0 that survive"""
    
    if d > l:           # vertex is outside collimator; no trajectories survive
        return 0.
    if d + 2*a0 <= l:   # vertex is sufficiently inside collimator that all trajectories sruvive
        return 1.
    
    c0 = 0
    if l+d < 2*a0:
        c0 = sqrt(1-((l+d)/(2.*a0))**2)  # cos phi below which no trajectories survive
        
    c1 = sqrt(1-((l-d)/(2.*a0))**2)      # cos phi above which all trajectories survive
    
    # integrate over partially-surviving portion
    return nintegrate(lambda c: plane_survival_fraction(l, d, sqrt(1-c**2)*a0), c0, c1) + 1-c1

def distrib_vertex_survival(l, a0):
    """Proportion of events surviving for vertices uniformly distributed over collimator"""
    return nintegrate((lambda d2: sphere_survival_fraction(l,sqrt(d2),a0)), 0, l**2)/l**2
    
def electron_Larmor(KE,B):
    """Electron Larmor radius [cm], for E in [keV] and B in [Gauss]"""
    p = sqrt(KE*KE+2*511.*KE);
    return 3.34*p/B

def plot_energy_survival():
    
    B = 364     # magnetic field [Gauss]
    l = 2.75    # electron collimator radius [cm]
    
    Emax = 1000
    
    g = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="energy [keV]", min=0, max=Emax),
        y=graph.axis.log(title="collimator pass fraction", min=0.001, max=1),
        key = graph.key.key(pos="tr"))
    
    gdat = [(KE, sphere_survival_fraction(l, 0, electron_Larmor(KE,B)), distrib_vertex_survival(l,electron_Larmor(KE,B)) ) for KE in unifrange(0, Emax, 500)]
    g.plot(graph.data.points(gdat,x=1,y=2,title="centered"),[graph.style.line([style.linewidth.Thick])])
    g.plot(graph.data.points(gdat,x=1,y=3,title="volume averaged"),[graph.style.line([style.linewidth.Thick,style.linestyle.dashed])])
    g.writetofile("Survival.pdf")
    
    
if __name__== "__main__":
    plot_energy_survival()
    
