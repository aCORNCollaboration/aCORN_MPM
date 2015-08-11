from math import *
from scipy import special
import numpy

from pyx import *
from pyx import style
from pyx import graph
from pyx.color import rgb, hsb
from pyx.graph.style import symbol

me = 511.    # electron mass, KeV/c^2
mp = 938000    # proton mass, KeV/c^2
Ee0 = 1293    # beta decay endpoint, keV

def proton_initial_momentum(Te):
    """Proton initial momenta in fast/slow groups as a function of electron KE"""
    Ee = Te + me
    pe = sqrt(Ee**2-me**2)    # electron momentum
    pn = Ee0 - Ee        # neutrino momentum
    return (pe + pn, pe - pn)
    
def proton_exit_velocity(Te, Em = 1.0):
    """Proton velocities at mirror exit after energy boost of Em for +/- groups"""
    pps = proton_initial_momentum(Te)
    return [ sqrt(2*(sqrt(p**2 + mp**2)-mp+Em)/mp) for p in pps]

def exit_velocity_ratio(Te, Em = 1.0):
    pvs = proton_exit_velocity(Te, Em)
    return pvs[0]/pvs[1]

if __name__ == "__main__":

    g = graph.graphxy(width=10,height=6,
            x=graph.axis.lin(title="$T_e$ [keV]", min = 0, max = Ee0-me),
            y=graph.axis.lin(title="$v_p^+/v_p^-$"))
    
    gdat = [ (Te, exit_velocity_ratio(Te)) for Te in numpy.linspace(0.1, Ee0-me-0.1, 400) ]
    g.plot(graph.data.points(gdat,x=1,y=2),[graph.style.line(lineattrs=[style.linewidth.thin,rgb.blue])])
        g.writetofile("exitv_ratio.pdf")

