#####################################################################
# Makefile for Michael P. Mendenhall's aCORN analysis code
#####################################################################


# assure correct shell is used
SHELL = /bin/sh
# apply implicit rules only for listed file types
.SUFFIXES:
.SUFFIXES: .c .cc .cpp .o
	 
# compiler command to use
CC = cc
CXX = g++

CXXFLAGS = -O3 --std=c++11 -fPIC `root-config --cflags` -I. -pedantic -Wall -Wextra \
	-IIOUtils -IROOTUtils -IBaseTypes -IMathUtils -ICalibration -IAnalysis -IStudies -IPhysics
LDFLAGS =  -L. -laCORN_MPM -lSpectrum -lMLP `root-config --libs` -lMathMore

ifdef PROFILER_COMPILE
	CXXFLAGS += -pg
	LDFLAGS += -pg
endif

ifdef UNBLINDED
	CXXFLAGS += -DUNBLINDED
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

VPATH = ./:IOUtils/:Physics/:ROOTUtils/:Standalone/:BaseTypes/:Analysis/

BaseTypes = TagCounter.o

IOUtils = ControlMenu.o OutputManager.o PathUtils.o QFile.o SMExcept.o strutils.o

Physics = BetaSpectrum.o ElectronBindingEnergy.o FloatErr.o NuclEvtGen.o

ROOTUtils = GraphUtils.o TChainScanner.o HistogramSequenceFitter.o

Analysis = RunSetScanner.o BaseDataScanner.o ReducedDataScanner.o SegmentSaver.o

objects = $(BaseTypes) $(IOUtils) $(Physics) $(ROOTUtils) $(Analysis)


all: libaCORN_MPM.a

libaCORN_MPM.a: $(objects)
	ar rs libaCORN_MPM.a $(objects)

# generic rule for everything else .cc linked against libaCORN_MPM
% : %.cc libaCORN_MPM.a
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

StandaloneObjs = ReducedToROOT PMT_Gainmatcher

standalone: $(StandaloneObjs)

#
# cleanup
#
.PHONY: clean
clean:
	-rm -f libaCORN_MPM.a
	-rm -f *.o
	-rm -rf *.dSYM

