#!/usr/bin/python

from PyxUtils import *
from numpy import *
import os
import cmath

class MapData:
	def __init__(self, fname, cols = {}):
		self.dat = array([[float(x) for x in l.split()] for l in open(fname,"r").readlines() if l[0] in "-+0123456789."])
		self.cols = cols
		self.cols["z"] = 0
		self.cols["t"] = 1
		self.cols["T"] = 2
		print "data array",self.dat.shape

def plot_processed(mapname, basedir):

	m = MapData(basedir+"/Processed_%s.txt"%mapname, {"Bx":3, "By": 4, "Bz": 5})

	g = graph.graphxy(width=10, height=8,
		x=graph.axis.lin(title="probe distance from top [cm]", min = 0, max = 250),
		y=graph.axis.lin(title="$B_z$ field [Gauss]", min=367, max=369))
		#key = graph.key.key(pos="tc",columns=2))

	g.plot(graph.data.points([(x[0], x[m.cols["Bz"]]) for x in m.dat], x=1, y=2, title=None),[graph.style.line()])
	g.writetofile("B_axial_%s.pdf"%mapname)

	gxy = graph.graphxy(width=10, height=8,
		x=graph.axis.lin(title="probe distance from top [cm]", min = 0, max = 250),
		y=graph.axis.lin(title="transverse field [milliGauss]", min=-80, max = 80),
		key = graph.key.key(pos="tc",columns=3))

	gxy.plot(graph.data.points([(x[0], 1000*x[m.cols["Bx"]]) for x in m.dat], x=1, y=2, title="$B_x$"), [graph.style.line([rgb.red])])
	gxy.plot(graph.data.points([(x[0], 1000*x[m.cols["By"]]) for x in m.dat], x=1, y=2, title="$B_y$"), [graph.style.line([rgb.blue])])
	gxy.plot(graph.data.points([(x[0], 1000*sqrt(x[m.cols["By"]]**2 + x[m.cols["Bx"]]**2)) for x in m.dat], x=1, y=2, title="$|B_T|$"), [graph.style.line([style.linestyle.dashed])])
	gxy.writetofile(basedir+"/B_xy_%s.pdf"%mapname)

	gProbes = graph.graphxy(width=10, height=8,
		x=graph.axis.lin(title="probe distance from top [cm]", min = 0, max = 250),
		#y=graph.axis.lin(title="probe tilt $\\sin \\theta$", min = -0.01, max = 0.01),
		y=graph.axis.lin(title="probe tilt $\\sin \\theta$", min = 0, max = 0.02),
		key = graph.key.key(pos="bc",columns=3))
	gProbes.plot(graph.data.points([(x[0], x[8]) for x in m.dat], x=1, y=2, title="$x$ Probe"), [graph.style.line([rgb.red])])
	gProbes.plot(graph.data.points([(x[0], x[11]) for x in m.dat], x=1, y=2, title="$y$ Probe"), [graph.style.line([rgb.blue])])
	if len(m.dat[0]) >= 15:
		gProbes.plot(graph.data.points([(x[0], sqrt(x[12]**2 + x[13]**2)) for x in m.dat], x=1, y=2, title="$z$ Probe"), [graph.style.line([style.linestyle.dashed])])
	gProbes.writetofile(basedir+"/ProbeTilt_%s.pdf"%mapname)

def plot_temperature(mapname, basedir):
	m = MapData(basedir+"/MapData_Bmap%s.txt"%mapname)
	
	g = graph.graphxy(width=10, height=8,
		x=graph.axis.lin(title="cart distance from top [cm]", min = 0, max = 270),
		y=graph.axis.lin(title="temperature [$^\\circ$C]"))

	g.plot(graph.data.points([(x[0], x[2]) for x in m.dat], x=1, y=2), [graph.style.line()])
	g.writetofile("Temperature_%s.pdf"%mapname)

def plot_raw_30(mapname, basedir, nskip = 1):
	m = MapData(basedir+"/MapData_Bmap%s.txt"%mapname)

	g = graph.graphxy(width=10, height=8,
		x=graph.axis.lin(title="probe distance from top [cm]", min = 0, max = 250),
		#y=graph.axis.lin(title="field reading [Gauss]"))
		y=graph.axis.lin(title="field reading [Gauss]", min = -5, max = 15))

	thcols = rainbow(13)
	for i in range(13)[::nskip]:
		g.plot(graph.data.points([(x[0], x[3 + 4*i + 0]) for x in m.dat], x=1, y=2), [graph.style.line([thcols[i]])])
		g.plot(graph.data.points([(x[0], x[3 + 4*i + 1]) for x in m.dat], x=1, y=2), [graph.style.line([style.linestyle.dashed,thcols[i]])])
		g.plot(graph.data.points([(x[0], x[3 + 4*i + 2]) for x in m.dat], x=1, y=2), [graph.style.line([style.linestyle.dotted,thcols[i]])])

	g.writetofile(basedir+"/RawScan_%s.pdf"%mapname)

	zcols = rainbow(m.dat.shape[0])
	iacols = [rgb.red, rgb.blue, rgb.green]
	gFFT = graph.graphxy(width=10, height=8,
			x=graph.axis.lin(title="probe distance from top [cm]", min = 0, max = 250),
			y=graph.axis.lin(title="Fourier components magnitude [milliGauss]", min = 0, max = 100),
			key = graph.key.key(pos="tc",columns=4))
	gPhase = graph.graphxy(width=10, height=8,
			x=graph.axis.lin(title="probe distance from top [cm]", min = 0, max = 250),
			y=graph.axis.lin(title="Fourier components phase [degrees]", min = 0, max = 360),
			key = graph.key.key(pos="tc",columns=4))

	for ia in range(3):
		axname = "xyz"[ia]
		if ia == 2 and nskip == 1:
			continue

		g = graph.graphxy(width=10, height=8,
			x=graph.axis.lin(title="probe angle [degrees]", min = 0, max = 360),
			y=graph.axis.lin(title="field reading [Gauss]", min = -0.1, max = 0.1))
	
		fftdat = []
		for iz in range(m.dat.shape[0]):
			# load theta scans
			scdat = [(30*i, m.dat[iz, 3 + 4*i + ia]) for i in range(13)[::nskip]]
			tdat = array([s[1] for s in scdat[:12/nskip]])
			z = m.dat[iz,0]			
			fftdat.append([z,] + list(fft.rfft(tdat)/len(tdat)))

			if not 25 < z < 225:
				continue

			# remove average values
			scdat = [(s[0], s[1]-fftdat[-1][1]) for s in scdat]
			g.plot(graph.data.points(scdat, x=1, y=2), [graph.style.line([zcols[iz]])])

		g.writetofile(basedir+"/RawScan_theta_%i_%s.pdf"%(ia,mapname))

		fncols = [[],[style.linestyle.dashed],[style.linestyle.dashdotted],[style.linestyle.dotted]]
		for i in range(4):
			if i+2 >= len(fftdat[0]):
				continue
			lsty = [graph.style.line(fncols[i] + [iacols[ia]])]
			gtitle = "$%i\\theta_%s$"%(i+1, axname) if i > 0 else "$\\theta_%s$"%axname
			# note factor of 2: real-FFT coefficient also has suppressed conjugate term.
			gFFT.plot(graph.data.points([(x[0],2*1000*abs(x[i+2])) for x in fftdat], x=1, y=2, title=gtitle), lsty)
			if i >= 2:
				continue
			dph = -90
			phdat = [(x[0], (cmath.phase(x[i+2])*180/pi + dph)%360, (cmath.phase(x[i+2])*180/pi - 90 + dph)%360, (cmath.phase(x[i+2])*180/pi + dph)%360) for x in fftdat]
			gPhase.plot(graph.data.points(phdat, x=1, y=2+ia, title=gtitle), lsty)

	gFFT.writetofile(basedir+"/RawScan_FFT_%s.pdf"%(mapname))
	gPhase.writetofile(basedir+"/RawScan_Phase_%s.pdf"%(mapname))

if __name__ == "__main__":
	basedir = "/home/mpmendenhall/Documents/aCORN/FieldMaps"
	
	# 1105: 90 degree, centered
	# 1107: 30 degree, off-center
	
	#plot_processed("1105")
	
	#plot_raw_30("1105", basedir, 3)
	plot_raw_30("1107", basedir)
	#plot_temperature("1105")
	#plot_temperature("1107")

