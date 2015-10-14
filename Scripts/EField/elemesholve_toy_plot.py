#!/usr/bin/python

from PyxUtils import *
from LinFitter import *
        
def sphere_toy_plot():
    lns = open("/home/mpmendenhall/Applications/elemesholve-bld/vtxdump.txt").readlines()
    gg = [[float(x) for x in l.split("\t")] for l in lns]
    gdat = [ [sqrt(x[0]**2 + x[1]**2 + x[2]**2), x[3]] for x in gg]
    gdat = [ x+[100*(x[1]*x[0] - 1.)] for x in gdat]
    
    g = graph.graphxy(width=10,height=6,
        x=graph.axis.log(title="radius [arb]"),
        y=graph.axis.log(title="potential $\\phi$ [arb]"))
    
    gR = graph.graphxy(width=10,height=6,
        x=graph.axis.log(title="radius [arb]"),
        y=graph.axis.lin(title="error $100 \\cdot (\\phi - \\phi_0)/\\phi_0$"))
     
    gR.plot(graph.data.points(gdat,x=1,y=3),[graph.style.symbol(symbol.circle, size = 0.04, symbolattrs=[style.linewidth.THin])])
    g.plot(graph.data.points(gdat,x=1,y=2),[graph.style.symbol(symbol.circle, size = 0.07)])
    g.plot(graph.data.function("y(x) = 1/x"),[graph.style.line(lineattrs=[rgb.red])])
    
    
    g.writetofile("twosphere.pdf")
    gR.writetofile("twosphere_resid.pdf")

def trace_compare():
    lns = open("../../Aux/Brian_v_Analytical_trace.txt").readlines()
    gdat = [ [float(x) for x in l.split()] for l in lns if l[0] != "#"]
    
    g = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="z [cm]", min=-6, max=6),
        y=graph.axis.lin(title="transverse field [V]", min = -1, max = 1.5),
        key=graph.key.key(pos="tl"))
    
    g.plot(graph.data.points(gdat, x=1, y=2, title = "Brian model"), [graph.style.line(lineattrs=[rgb.red])])
    g.plot(graph.data.points(gdat, x=1, y=5, title = "Analytical, r=50um"), [graph.style.line(lineattrs=[rgb.blue])])
    g.writetofile("FieldTrace.pdf")
    
def integral_compare():
    lns = open("../../Aux/elemesholve_analytical_transverse_integrated_fields.txt").readlines()
    gdat = [ [float(x) for x in l.split()] for l in lns if l[0] != "#"]
    
    lnsM = open("../../Aux/elemesholve_transverse_integrated_fields.txt").readlines()
    gdatM = [ [float(x) for x in l.split()] for l in lnsM if l[0] != "#"]
    
    lnsB = open("../../Aux/Brian_transverse_integrated_fields.txt").readlines()
    gdatB = [ [float(x) for x in l.split()] for l in lnsB if l[0] != "#"]
    
    g = graph.graphxy(width=10,height=6,
        x=graph.axis.lin(title="radius [cm]", min=0, max=4),
        y=graph.axis.lin(title="integrated transverse field $\\int E_r dz$ [V]", min=0, max=8),
        key=graph.key.key(pos="tl"))
    
    g.plot(graph.data.points(gdatM, x=1, y=2, title = "Mesh inside mirror"), [graph.style.symbol(symbol.circle, size = 0.05, symbolattrs=[rgb.green])])
    #g.plot(graph.data.points(gdatM, x=1, y=3, title = "Mesh outside mirror"), [graph.style.symbol(symbol.triangle, size = 0.05, symbolattrs=[rgb.green])])
    
    inattrs = [style.linestyle.dashed, style.linewidth.thick]
    g.plot(graph.data.points(gdatB, x=1, y=5, title = "Analytical inside mirror"), [graph.style.line(lineattrs=inattrs+[rgb.blue])])
    g.plot(graph.data.points(gdatB, x=1, y=2, title = "Brian inside mirror"), [graph.style.line(lineattrs=inattrs+[rgb.red])])
    
    oattrs = [style.linewidth.thick]
    g.plot(graph.data.points(gdatB, x=1, y=6, title = "Analytical outside mirror"),[graph.style.line(lineattrs=oattrs+[rgb.blue])])
    g.plot(graph.data.points(gdatB, x=1, y=3, title = "Brian outside mirror"), [graph.style.line(lineattrs=oattrs+[rgb.red])])
    
    #g.plot(graph.data.points(gdatB, x=1, y=7, title = "Analytical total"), [graph.style.line(lineattrs=[rgb.blue])])
    g.plot(graph.data.points(gdatB, x=1, y=4, title = "Brian total"), [graph.style.line(lineattrs=[rgb.red])])
    
    LF = LinearFitter([polyterm(n) for n in range(6)][1:])
    LF.fit([p for p in gdatB if 0 <= p[0] <= 3.95], cols = (0,3))
    g.plot(graph.data.points(LF.fittedpoints(), x=1, y=3, title = None), [graph.style.line()])
    print LF.toLatex()
    print LF.coeffs
    
    g.writetofile("Etransverse.pdf")
    
if __name__=="__main__":
    trace_compare()
    integral_compare()
    