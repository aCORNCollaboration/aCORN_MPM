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

# check for necessary environment variables
ifndef MPMUTILS
$(error MPMUTILS is not set)
endif

# apply implicit rules only for listed file types
.SUFFIXES:
.SUFFIXES: .c .cc .cpp .o

# compiler flags
CXXFLAGS = -O3 --std=c++11 -fPIC `root-config --cflags` -I. -pedantic -Wall -Wextra \
	-I${MPMUTILS}/GeneralUtils/ -I${MPMUTILS}/Physics/ -I${MPMUTILS}/Matrix/ -I${MPMUTILS}/ROOTUtils/ \
	-IIOUtils -IROOTUtils -IBaseTypes -ICalibration -IAnalysis -IStudies -IPhysics
LDFLAGS +=  -L. -L${MPMUTILS}/GeneralUtils/ -L${MPMUTILS}/ROOTUtils/ \
	-laCORN_MPM -lMPMROOTUtils -lMPMGeneralUtils -lSpectrum -lMLP `root-config --libs` -lMathMore -lgsl -lblas

ifdef PROFILER_COMPILE
	CXXFLAGS += -pg
	LDFLAGS += -pg
endif

#
# things to build
#

VPATH = ./:IOUtils/:Physics/:Standalone/:BaseTypes/:Analysis/:${MPMUTILS}/GeneralUtils/:${MPMUTILS}/Physics/

IOUtils = sqlite3.o TextTableReader.o

Physics = PolarizedBetaAsym.o ElectronBindingEnergy.o NuclEvtGen.o \
	UnpolarizedBeta.o UnpolarizedNeutronDecay.o Collimator.o

ROOTUtils = SQLite_Helper.o

Analysis = AcornCalibrator.o AcornDB.o BaseDataScanner.o Positioner.o PMTsPlugin.o PluginInterpolator.o \
	ReducedDataScanner.o ResponseMatrix.o RunAccumulator.o RunSetScanner.o SimWishboneSmearer.o \
	SourceCalPlugin.o WishboneFit.o WishbonePlugin.o

objects = $(IOUtils) $(Physics) $(ROOTUtils) $(Analysis) aCornCompileVersion.o
StandaloneObjs = BetaSpectrometerScanner ReducedToROOT PMT_Gainmatcher WishboneScanner SpectrumShape

all: libaCORN_MPM.a $(StandaloneObjs)

libaCORN_MPM.a: $(objects)
	ar rs libaCORN_MPM.a $(objects)

.PHONY: aCornCompileVersion.o
aCornCompileVersion.o: 
	$(CXX) -c $(CXXFLAGS) -DGIT_SHA=$(shell git rev-parse -q HEAD) BaseTypes/aCornCompileVersion.cc -o aCornCompileVersion.o

# generic rule for everything else .cc linked against libaCORN_MPM
% : %.cc libaCORN_MPM.a
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

# documentation
.PHONY: doc
doc:
	doxygen Doxyfile

# cleanup
.PHONY: clean
clean:
	-rm -f *.a *.o *.dSYM
	-rm -f $(StandaloneObjs)
