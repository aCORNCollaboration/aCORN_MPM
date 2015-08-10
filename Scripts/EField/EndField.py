from math import *
from scipy import special
import numpy

from pyx import *
from pyx import style
from pyx import graph
from pyx.color import rgb, hsb
from pyx.graph.style import symbol

class BesselCalcs:
	def __init__(self):
		self.j0n = [0,]+list(special.jn_zeros(0, 1000))
		self.cns = [0,] + [self.cn(n+1) for n in range(100)]
	def cn(self,n):
		return 2/(self.j0n[n]*special.jv(1,self.j0n[n]))
	def V(self, r, z):
		s = 0
		for (n,c) in enumerate(self.cns):
			if n==0:
				continue
			s += self.cns[n]*exp(-self.j0n[n]*z)*special.jv(0,self.j0n[n]*r)
		return s

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

