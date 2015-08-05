#!/usr/bin/python

from ElectronSurvival import *
from scipy import special
from scipy import integrate
import numpy

me = 511.       # electron mass, KeV/c^2
mp = 938000     # proton mass, KeV/c^2
Ee0 = 1293      # beta decay endpoint, keV

def proton_initial_momentum(Te):
        """Proton initial momenta in fast/slow groups as a function of electron KE"""
        Ee = Te + me
        pe = sqrt(Ee**2-me**2)  # electron momentum
        pn = Ee0 - Ee           # neutrino momentum
        return (pe + pn, pe - pn)
        
def proton_exit_velocity(Te, Em = 1.0):
        """Proton velocities at mirror exit after energy boost of Em for +/- groups"""
        pps = proton_initial_momentum(Te)
        return [ sqrt(2*(sqrt(p**2 + mp**2)-mp+Em)/mp) for p in pps]

def exit_velocity_ratio(Te, Em = 1.0):
        pvs = proton_exit_velocity(Te, Em)
        return pvs[0]/pvs[1]
    
    
    
# l: collimator diameter
# d: distance of origin point from collimator
# a: Larmor radius
    
def proton_Larmor(KEe,B):
    """Proton max Larmor radius [cm], given electron KE in [keV] and B in [Gauss]"""
    Enu = Ee0 - me - KEe   # neutrino energy = neutrino momentum * c
    pp = Enu            # max proton transverse momentum = neutrino momentum
    return 3.34*pp/B    # Larmor radius [cm]
    
def shifted_survival_fraction(l, d0, a, da):
    """Survival fraction after momentum shift from initial center d0=(x0,y0) & larmor radius a, by da=(dx,dy)"""
    
    d1 = ( d0[0] + da[0], d0[1] + da[1] )       # center of new orbits
    r1 = sqrt( d1[0]**2 + d1[1]**2 )            # radius of new orbits center
    if r1 > l:           # new orbit centers outside collimator: none pass
        return 0

    mda2 = da[0]**2 + da[1]**2                  # magnitude^2 of momentum change
    if not mda2:         # zero momentum change: hard cut on whether single orbit passes
        return r1 + a < l
    
    mda = sqrt(mda2)                            # magnitude of momentum change
    
    if not a:           # zero initial momentum: hard cut on whether single orbit passes
        return r1 + mda < l
    
    # minimum cos theta on original orbit that will pass collimator
    cthmin = (a**2 + mda2 - (l-r1)**2)/(2*a*mda)
    if cthmin >= 1:
        return 0
    if cthmin <= -1:
        return 1
    return acos(cthmin)/pi
    
    
    
def shifted_distrib_vertex_survival(l, r0, a0, da):
    """Survival fraction for uniform initial spiral centers in a ring at r0, initial radius a0"""
    if r0 + a0 + 2*da < l:
        return 1
    I = integrate.quad(lambda th: shifted_survival_fraction(l, (r0*cos(th), r0*sin(th)), a0, (da, 0)), 0, pi, epsrel = 1e-2)
    return I[0]/pi




def shifted_sphere_c_range(l, r0, a0, da):
    csqmin = 0
    csqmid = 0.5
    csqmax = 1
    # r0 - da - |da - a*| < l ?
    if 0 < a0:
        csqmid = 1-((l-r0)/a0)**2
        csqmin = 1-((l-r0+2*da)/a0)**2
        if csqmin > csqmid:
            csqmin = csqmid
        csqmax = 1-((l-r0-2*da)/a0)**2
        if (l-r0-2*da) < 0:
            csqmax = 1
    if csqmax < 0:
        csqmax = 0
    if csqmin < 0:
        csqmin = 0
    if csqmid < 0:
        csqmid = 0
    
    if not csqmin <= csqmid <= csqmax:
        print l, r0, a0, da, csqmin, csqmid, csqmax
    assert csqmin <= csqmid <= csqmax
    
    return (sqrt(csqmin), sqrt(csqmid), sqrt(csqmax))

def shifted_sphere_theta_min(l, r0, a0, da):
    thmin = 0
    if r0 and da:
        u = ((l-abs(a0 - da))**2-da**2-r0**2)/(2*r0*da)
        if -1 < u < 1:
            thmin = acos(u)
    return thmin

def shifted_sphere_theta_max(l, r0, a0, da):
    thmax = pi
    if r0 and da:
        u = ((l-(a0 + da))**2-da**2-r0**2)/(2*r0*da)
        if -1 < u < 1:
            thmax = acos(u)
    return thmax

def shifted_sphere_vertex_survival(l, r0, a0, da):
    """Shifted survival fraction for points originating on Larmor spirals at a common center r0 from collimator center, over a sphere with max Larmor radius a0"""
    
    crange = shifted_sphere_c_range(l, r0, a0, da)
    if crange[2] == 0:
        return 1
    
    f = 0
    
    if 1:
        Ia = integrate.dblquad(lambda th, c: shifted_survival_fraction(l, (r0*cos(th), r0*sin(th)), a0*sqrt(1-c*c), (da,0)),
                            crange[0], crange[1], 
                            (lambda c: shifted_sphere_theta_min(l, r0, a0*sqrt(1-c*c), da)),
                            (lambda c: shifted_sphere_theta_max(l, r0, a0*sqrt(1-c*c), da)),
                            epsrel = 1e-3, epsabs = 1e-4)
        Ia1 = integrate.quad(lambda c: pi-shifted_sphere_theta_max(l, r0, a0*sqrt(1-c*c), da), crange[0], crange[1], epsrel = 1e-3, epsabs = 1e-4)
        
        Ib = integrate.dblquad(lambda th, c: shifted_survival_fraction(l, (r0*cos(th), r0*sin(th)), a0*sqrt(1-c*c), (da,0)),
                            crange[1], crange[2], 
                            (lambda c: shifted_sphere_theta_min(l, r0, a0*sqrt(1-c*c), da)),
                            (lambda c: shifted_sphere_theta_max(l, r0, a0*sqrt(1-c*c), da)),
                            epsrel = 1e-3, epsabs = 1e-4)
        Ib1 = integrate.quad(lambda c: pi-shifted_sphere_theta_max(l, r0, a0*sqrt(1-c*c), da), crange[1], crange[2], epsrel = 1e-3, epsabs = 1e-4)
        
        f = (Ia[0]+Ia1[0]+Ib[0]+Ib1[0])/pi + (1-crange[2])
    else:
        Ia = integrate.dblquad(lambda th, c: shifted_survival_fraction(l, (r0*cos(th), r0*sin(th)), a0*sqrt(1-c*c), (da,0)),
                            crange[0], crange[2], 
                            (lambda c: shifted_sphere_theta_min(l, r0, a0*sqrt(1-c*c), da)),
                            (lambda c: shifted_sphere_theta_max(l, r0, a0*sqrt(1-c*c), da)),
                            epsrel = 1e-3, epsabs = 1e-4)
        Ia1 = integrate.quad(lambda c: pi-shifted_sphere_theta_max(l, r0, a0*sqrt(1-c*c), da), crange[0], crange[2], epsrel = 1e-3, epsabs = 1e-4)
        f = (Ia[0]+Ia1[0])/pi + (1-crange[2])
    
    assert 0 <= f <= 1.0001
    return f
    
def shifted_sphere_vertex_survival_plots(l, r0, a0, da):
    
    crange = shifted_sphere_c_range(l, r0, a0, da)
    print crange
    
    g = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="$\\theta$", min=0, max=pi),
        y=graph.axis.lin(title="proton survival fraction", min=0.0, max=1),
        key = graph.key.key(pos="tl"))
    gc = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="$c$", min=crange[0], max=crange[2]),
        y=graph.axis.lin(title="proton survival fraction"),
        key = graph.key.key(pos="tl"))

    rsty = [style.linestyle.solid, style.linestyle.dashed, style.linestyle.dotted, style.linestyle.dashdotted]
    gdatc = []
    for (n,c) in enumerate(list(unifrange(crange[0],crange[1],50)) + list(unifrange(crange[1],crange[2],50))):
        thmn = shifted_sphere_theta_min(l, r0, a0*sqrt(1-c*c), da)
        thmx = shifted_sphere_theta_max(l, r0, a0*sqrt(1-c*c), da)
        gdat = [(th, shifted_survival_fraction(l, (r0*cos(th), r0*sin(th)), a0*sqrt(1-c*c), (da,0))) for th in unifrange(thmn,thmx,200)]
        gdatc.append((c, sum([x[1] for x in gdat])*(thmx-thmn)/pi/len(gdat) + (pi - thmx)/pi))
        g.plot(graph.data.points(gdat,x=1,y=2,title=None),[graph.style.line([style.linewidth.THIn, rsty[n%4]])])
    
    g.writetofile("SphereIntegral.pdf")   
    
    gc.plot(graph.data.points(gdatc,x=1,y=2,title=None),[graph.style.line([style.linewidth.thick])])
    gc.writetofile("SphereIntegral_c.pdf")
    
    
 
def uvtx_splitpts(l, a0, da): 
    #splitpts = [0, da, l-(a0 + 2*da), l-2*da, l-da, l-a0, l-(a0+da), l-(a0-da), l-(a0-2*da), l+a0, l+da, l+a0+da, l+a0+2*da, 1.5*l]
    splitpts = [0, da, a0, abs(l-(a0 + 2*da)), abs(l-a0), abs(l-da), abs(l-abs(a0-da)), l, 1.5*l]
    splitpts.sort()
    splitpts = [ p for p in splitpts if 0 <= p <= 1.5*l]
    return splitpts

def uniform_vertex_survival(l, a0, da):
    """Total survival fraction for uniformly distributed initial vertices (normalized to collimation in vertex plane)"""
    
    splitpts = uvtx_splitpts(l, a0, da) 
    print (l,a0,da), splitpts
    
    I = 0
    for n in range(len(splitpts)-1):
        l0 = 1.0001*splitpts[n]+0.00001
        l1 = 0.9999*splitpts[n+1]-0.00001
        if l1 <= l0:
            continue
        #return nintegrate((lambda d2: shifted_sphere_vertex_survival(l, sqrt(d2), a0, da)), 0, 2*l**2)/l**2
        #return integrate.quad(lambda d2: shifted_sphere_vertex_survival(l, sqrt(d2), a0, da), 0, 2*l**2, epsrel = 1e-2, epsabs = 1e-3)[0]/l**2        
        I += integrate.quad(lambda d: d*shifted_sphere_vertex_survival(l, d, a0, da), l0, l1, epsrel = 1e-2, epsabs = 1e-3)[0]/(2*l)

    return I

def uniform_vertex_survival_plot(l, a0, da):
    
    llmax = 1.5*l**2
    
    g = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="$l$", min=0, max = sqrt(llmax/l**2)),
        y=graph.axis.lin(title="proton survival fraction", min = 0))
   
   
    splitpts = uvtx_splitpts(l, a0, da)
    print (l,a0,da), splitpts
    rsty = [style.linestyle.solid, style.linestyle.dashed, style.linestyle.dotted, style.linestyle.dashdotted]
    
    for n in range(len(splitpts)-1):
        l0 = 1.0001*splitpts[n]
        l1 = 0.9999*splitpts[n+1]
        if l1 <= l0:
            continue
        gdat = [ (ll/l, l*shifted_sphere_vertex_survival(l, ll, a0, da)) for ll in unifrange(l0,l1,50)]
        g.plot(graph.data.points(gdat,x=1,y=2),[graph.style.line([style.linewidth.thick, rsty[n%4]])])

    g.writetofile("UniformVertex.pdf")   
    
    
    
    
    
    
    
def plot_proton_survival():
    B = 364     # magnetic field [Gauss]
    l = 4.0     # proton collimator radius [cm]
    
    Emax = Ee0-me
    
    g = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="electron kinetic energy [keV]", min=0, max=Emax),
        y=graph.axis.log(title="proton survival fraction", min=0.02, max=1),
        key = graph.key.key(pos="tl"))
    
    rsty = [style.linestyle.solid, style.linestyle.dashed, style.linestyle.dotted, style.linestyle.dashdotted]
    
    for (n,rp) in enumerate([0,1,2,3]):        # initial proton radius
        gdat = [(KE, sphere_survival_fraction(l, rp, proton_Larmor(KE,B))) for KE in unifrange(0,Emax,100) ]
        g.plot(graph.data.points(gdat,x=1,y=2,title="$r_0 = %g$ cm"%rp),[graph.style.line([style.linewidth.thick, rsty[n]])])
    g.writetofile("ProtonSurvival.pdf")
    
def plot_shifted_survival():
    B = 364     # magnetic field [Gauss]
    l = 4.0     # proton collimator radius [cm]
    Emax = Ee0-me
    
    g = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="electron kinetic energy [keV]", min=0, max=Emax),
        y=graph.axis.log(title="proton survival fraction", min=0.01, max=2),
        key = graph.key.key(pos="br"))
    
    rsty = [style.linestyle.solid, style.linestyle.dashed, style.linestyle.dotted, style.linestyle.dashdotted]
    adp = 3.34*30/B # momentum Larmor offset (cm) for 30 keV/c = 1 eV us / cm
    #print "a(dp) =",adp
    
    #shifted_sphere_vertex_survival_plots(l, 1.6, proton_Larmor(450,B), adp)
    #shifted_sphere_vertex_survival_plots(4.0, 3.46410161514, 7.17549450549, 0.594729125864)
    #uniform_vertex_survival_plot(l, proton_Larmor(200,B), 0.2)
    #return

    if 0:
        for (n,rp) in enumerate([0,1,2,3]):        # initial proton radius
            print rp
            #gdat = [(KE, shifted_distrib_vertex_survival(l, rp, proton_Larmor(KE,B), adp)) for KE in unifrange(0, Emax, 500) ]
            gdat = [(KE, shifted_sphere_vertex_survival(l, rp, proton_Larmor(KE,B), adp)) for KE in unifrange(0, Emax, 30) ]
            g.plot(graph.data.points(gdat,x=1,y=2,title="$r_0 = %g$ cm"%rp),[graph.style.line([style.linewidth.thick, rsty[n%4]])])
        #for (n,rp) in enumerate([0,1,2,3]):        # initial proton radius
        #    gdat = [(KE, shifted_sphere_vertex_survival(l, rp, proton_Larmor(KE,B), exit_velocity_ratio(KE)*adp)) for KE in unifrange(0,Emax,100) ]
        #    g.plot(graph.data.points(gdat,x=1,y=2,title=None),[graph.style.line([style.linewidth.thick, rsty[n]])])
    
    if 1:
        gdat = []
        for KE in unifrange(0,Emax,15):
            vexit = proton_exit_velocity(KE)
            da0 =  3e-4/vexit[0]
            da1 =  3e-4/vexit[1]
            a0 = proton_Larmor(KE,B)
            acc0 = uniform_vertex_survival(l, a0, da0)
            acc1 = uniform_vertex_survival(l, a0, da1)
            gdat.append((KE, acc0, acc1, (acc0-acc1)/(acc0 + acc1)))
            print gdat[-1], da0, da1, vexit
        g.plot(graph.data.points(gdat,x=1,y=2,title=None),[graph.style.line([style.linewidth.thick])])
        g.plot(graph.data.points(gdat,x=1,y=3,title=None),[graph.style.line([style.linewidth.THin])])
    
    g.writetofile("ShiftedSurvival.pdf")
    
    if 1:
        gA = graph.graphxy(width=10,height=6,
            x=graph.axis.lin(title="electron kinetic energy [keV]", min=0, max=Emax),
            y=graph.axis.lin(title="false asymmetry"))
        gA.plot(graph.data.points(gdat,x=1,y=4,title=None),[graph.style.line([style.linewidth.thick])])
        gA.writetofile("FalseAsym.pdf")
     
if __name__== "__main__":
    #plot_proton_survival()
    plot_shifted_survival()
