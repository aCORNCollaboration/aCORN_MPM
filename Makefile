#####################################################################
# Makefile for Michael P. Mendenhall's aCORN analysis code
#####################################################################
# This file was produced under the employ of the United States Government,
# and is consequently in the PUBLIC DOMAIN, free from all provisions of
# US Copyright Law (per USC Title 17, Section 105).
#
# -- Michael P. Mendenhall, 2015


# assure correct shell is used
SHELL = /bin/sh
# apply implicit rules only for listed file types
.SUFFIXES:
.SUFFIXES: .c .cc .cpp .o
	 
# compiler flags
CXXFLAGS = -O3 --std=c++11 -fPIC `root-config --cflags` -I. -pedantic -Wall -Wextra \
	-I${MPMUTILS}/GeneralUtils/ -I${MPMUTILS}/Matrix/ -I${MPMUTILS}/ROOTUtils/ \
	-IIOUtils -IROOTUtils -IBaseTypes -ICalibration -IAnalysis -IStudies -IPhysics
LDFLAGS =  -L. -L${MPMUTILS}/GeneralUtils/ -L${MPMUTILS}/ROOTUtils/ \
	-laCORN_MPM -lMPMROOTUtils -lMPMGeneralUtils -lSpectrum -lMLP `root-config --libs` -lMathMore -lgsl -lblas

ifdef PROFILER_COMPILE
	CXXFLAGS += -pg
	LDFLAGS += -pg
endif

ifdef PUBLICATION_PLOTS
	CXXFLAGS += -DPUBLICATION_PLOTS
endif

ifdef TSPECTRUM_USES_DOUBLE
	CXXFLAGS += -DTSPECTRUM_USES_DOUBLE
endif

#
# things to build
#

VPATH = ./:IOUtils/:Physics/:Standalone/:BaseTypes/:Analysis/:${MPMUTILS}/GeneralUtils/

IOUtils = sqlite3.o

Physics = BetaSpectrum.o ElectronBindingEnergy.o NuclEvtGen.o UnpolarizedBeta.o aSpectrum.o Collimator.o

ROOTUtils = SQLite_Helper.o

Analysis = AcornCalibrator.o AcornDB.o BaseDataScanner.o Positioner.o PMTsPlugin.o PluginInterpolator.o \
	ReducedDataScanner.o RunAccumulator.o RunSetScanner.o SourceCalPlugin.o WishboneFit.o WishbonePlugin.o

objects = $(IOUtils) $(Physics) $(ROOTUtils) $(Analysis)


all: libaCORN_MPM.a

libaCORN_MPM.a: $(objects)
	ar rs libaCORN_MPM.a $(objects)

# generic rule for everything else .cc linked against libaCORN_MPM
% : %.cc libaCORN_MPM.a
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

StandaloneObjs = BetaSpectrometerScanner ReducedToROOT PMT_Gainmatcher WishboneScanner

standalone: $(StandaloneObjs)

#
# cleanup
#
.PHONY: clean
clean:
	-rm -f libaCORN_MPM.a
	-rm -f *.o
	-rm -rf *.dSYM
	-rm -f ReducedToROOT  PMT_Gainmatcher SpectrumShape  WishboneScanner CalibSpectra


