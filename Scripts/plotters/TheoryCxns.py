#!/usr/bin/python

from math import *
import numpy
from pyx import *
from pyx import style
from pyx import graph
from pyx.color import rgb, hsb
from pyx.graph.style import symbol
from LinFitter import *

def rainbow(n, b=1.0):
    return [ hsb((1.0*x)/n,1,b) for x in range(n) ]
    
def rainbowDict(keys, b=1.0):
    n = len(keys)
    knew = [k for k in keys]
    knew.sort()
    return dict([ (k,hsb((1.0*x)/n,1,b)) for (x,k) in enumerate(knew) ])
        
l = abs(-1.2723)    # +/-0.0023, |lambda| PDG 2014, 
mu = 2.792847356-(-1.91304273)
m = 511.00          # keV, m_e
M = 939565.379      # keV, +/- 21eV, neutron mass, PDG 2014 
D = 1293.35         # keV, M_p-M_n
alpha = 1./137.

a0 = (1-l**2)/(1+3*l**2)

def beta(KE):
    if KE <= 0:
        return 1
    return sqrt(KE*KE+2*m*KE)/(m+KE)

def a_ideal(E,theta):
    return 1 + a0*beta(E-m)*cos(theta)
    
def bilenkii_1959_a_cxn(E,theta):
    a = a0 + ( 4*l*(1+l**2)*(l+mu)*D/M
                + (1-l**2)*(1+l**2+2*l*mu)*m**2/(M*E)
                -(8*l*(1+l**2)*mu + 3*(1+6*l**2+9*l**4))*E/M ) / (1+3*l**2)**2
    b = 3*(l**2-1)/(1+3*l**2)*E/M
    x = beta(E-m)*cos(theta)
    return 1 + a*x + b*x**2

def hmg(E):
    b = beta(E-m)
    a = atanh(b)
    return 4*(a/b-1)*(1-b**2+(D-E)/(8*E))*(D-E)/(3*E*b**2) + a/b*(2-2*b**2-(D-E)**2/(6*E**2))
    
def plot_bilenkii_a_cxn():
    npts = 100
    gdat = []
    for i in range(npts+1):
        E = m + i*(D-m)/npts
        a0p = a_ideal(E,0)
        a0m = a_ideal(E,pi)
        ap = bilenkii_1959_a_cxn(E,0)
        am = bilenkii_1959_a_cxn(E,pi)
        
        a0 = (a0p-a0m)/(a0p+a0m)
        a = (ap-am)/(ap+am)
        gdat.append([E-m,(a-a0)/a0*100,])

    g = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="kinetic energy [keV]",min=0,max=D-m),
        y=graph.axis.lin(title="$(a_{\\rm meas}-a_0)/a_0$ [\\%]"))
        
    g.plot(graph.data.points(gdat,x=1,y=2),[graph.style.line([style.linewidth.THIck])])
    g.writetofile("Bilenkii_1959_a_cxn.pdf")

def gluck_1993_r_enu(E,theta):
    c = cos(theta)
    x = (E-m)/(D-m)
    r_enu = (0.002 + 0.014*x) + (-0.009 + 0.179*x)*c
    print "r_enu",x,c,r_enu
    return r_enu

def gluck_1993_r_e(E):
    x = (E-m)/(D-m)
    xr_e = [    0.1,    0.2,    0.3,    0.4,    0.5,    0.6,    0.7,    0.8,    0.9,    0.95 ]
    r_e = [     1.82,   1.74,   1.65,   1.53,   1.40,   1.25,   1.07,   0.84,   0.50,   0.21 ]
    return numpy.interp(x, xr_e, r_e)
    
def plot_gluck1993_a_cxn():
    npts = 100
    gdat = []
    for i in range(npts+1)[1:]:
        E = m + i*(D-m)/npts
        
        a0p = a_ideal(E,0)
        a0m = a_ideal(E,pi)
        r_e = 0.01*gluck_1993_r_e(E)
        ap = a0p*(1 + r_e + 0.01*gluck_1993_r_enu(E,0))
        am = a0m*(1 + r_e + 0.01*gluck_1993_r_enu(E,pi))
        
        print E, a0p, a0m
        print "\t", ap, am
        print "\t", gluck_1993_r_enu(E,0), gluck_1993_r_enu(E,pi)
        
        a0 = (a0p-a0m)/(a0p+a0m)
        a = (ap-am)/(ap+am)
        gdat.append([E-m,(a-a0)/a0*100, 100*alpha*hmg(E)/(2*pi)])
        
    g = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="kinetic energy [keV]",min=0,max=D-m),
        y=graph.axis.lin(title="$(a_{\\rm meas}-a_0)/a_0$ [\\%]"),
        key=graph.key.key(pos="bl"))
        
    g.plot(graph.data.points(gdat,x=1,y=2,title='Gl\\"uck'),[graph.style.line([style.linewidth.THick])])
    g.plot(graph.data.points(gdat,x=1,y=3,title="Garc\\'ia-Maya"),[graph.style.line([style.linewidth.Thick, style.linestyle.dotted])])
    g.writetofile("Gluck_1993_a_cxn.pdf")
    #$\\frac{\\alpha}{2\\pi}(h-g)$"
    
def gluck_1993_cnxs():
    xr_e = [    0.1,    0.2,    0.3,    0.4,    0.5,    0.6,    0.7,    0.8,    0.9,    0.95 ]
    r_e = [     1.82,   1.74,   1.65,   1.53,   1.40,   1.25,   1.07,   0.84,   0.50,   0.21 ]
    
    # Table V, [[c,   [r_enu(x)...]], ... ]
    r_enu = [   [0.9, [ 0.03,   0.04,   0.06,   0.08,   0.10,   0.12,   0.14,   0.16 ]],
                [0.7, [ 0.03,   0.04,   0.05,   0.07,   0.08,   0.10,   0.11,   0.13 ]],
                [0.5, [ 0.02,   0.03,   0.04,   0.05,   0.06,   0.07,   0.08,   0.09 ]],
                [0.3, [ 0.02,   0.02,   0.03,   0.03,   0.04,   0.04,   0.05,   0.06 ]],
                [0.1, [ 0.01,   0.01,   0.01,   0.02,   0.02,   0.02,   0.02,   0.03 ]],
                [-0.1,[ 0.01,   0.00,   0.00,   -0.00,  -0.00,  -0.00,  -0.01,  -0.00 ]],
                [-0.3,[ -0.00,  -0.01,  -0.01,  -0.02,  -0.02,  -0.03,  -0.03,  -0.03 ]],
                [-0.5,[ -0.01,  -0.02,  -0.02,  -0.05,  -0.04,  -0.05,  -0.05,  -0.06 ]],
                [-0.7,[ -0.01,  -0.02,  -0.02,  -0.05,  -0.06,  -0.07,  -0.08,  -0.08 ]],
                [-0.9,[ -0.02,  -0.03,  -0.05,  -0.06,  -0.08,  -0.09,  -0.10,  -0.11 ]] ]
                
    # w = F(E)*[1+a*beta*c]*[1 + 0.01*r_e(x) + 0.01*r_enu(x,c)]
    
    gx = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title='$x$',min=0,max=1),
        y=graph.axis.lin(title='Gl\\"uck $r_e(x) + r_{e\\nu}(x,c)$'))
    
    gc = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title='$c$',min=-1,max=1),
        y=graph.axis.lin(title='Gl\\"uck $r_{e\\nu}(x,c)$'),
        key=graph.key.key(pos="tl",columns=2))
    
    ccols = rainbowDict([r[0] for r in r_enu])
    xcols = rainbowDict(xr_e[1:-1])
    
    cdat = [[x,[]] for x in xr_e[1:-1]]
    for (c, r) in r_enu:
        for (n,rr) in enumerate(r):
            cdat[n][1].append([c,rr])
        gdat = [(xr_e[n+1], r_e[n+1]+rr) for (n,rr) in enumerate(r)]
        gx.plot(graph.data.points(gdat,x=1,y=2),[graph.style.symbol(symbol.circle,size=0.1,symbolattrs=[ccols[c]]), graph.style.line(lineattrs=[ccols[c]])])
    
    gx.writetofile("Gluck_1993_a_cxn_x.pdf")
    
    paramdat = []
    for (x,gdat) in cdat:
        LF = LinearFitter(terms=[(lambda z: 1), (lambda z: z)])
        LF.fit(gdat)
        paramdat.append([x,LF.coeffs[0],LF.coeffs[1]])
        gc.plot(graph.data.points(gdat,x=1,y=2,title="$x = %g$"%x),[graph.style.symbol(symbol.circle,size=0.1,symbolattrs=[xcols[x]])])
        gc.plot(graph.data.points(LF.fitcurve(-1,1),x=1,y=2,title=None),[graph.style.line(lineattrs=[xcols[x]])])
    gc.writetofile("Gluck_1993_a_cxn_c.pdf")
    
    gcparam = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title='$x$',min=0,max=1),
        y=graph.axis.lin(title='$r_{e\\nu}$ parametrization term',min=0, max=0.2),
        key=graph.key.key(pos="tl"))
        
    LF.fit(paramdat,cols=(0,1))
    gcparam.plot(graph.data.points(paramdat,x=1,y=2,title="$r_{e\\nu}^0(x) \\approx %.3f+%.3f \\cdot x$"%tuple(LF.coeffs)),[graph.style.symbol(symbol.circle,size=0.1,symbolattrs=[rgb.red])])
    gcparam.plot(graph.data.points(LF.fitcurve(0,1),x=1,y=2,title=None),[graph.style.line(lineattrs=[rgb.red])])
    
    LF.fit(paramdat,cols=(0,2))
    gcparam.plot(graph.data.points(paramdat,x=1,y=3,title="$r_{e\\nu}^1(x) \\approx %.3f+%.3f \\cdot x$"%tuple(LF.coeffs)),[graph.style.symbol(symbol.circle,size=0.1,symbolattrs=[rgb.blue])])
    gcparam.plot(graph.data.points(LF.fitcurve(0,1),x=1,y=2,title=None),[graph.style.line(lineattrs=[rgb.blue])])
    
    gcparam.writetofile("Gluck_1993_a_cxn_cparam.pdf")
    
if __name__ == "__main__":
    #plot_bilenkii_a_cxn()
    plot_gluck1993_a_cxn()
    #gluck_1993_cnxs()