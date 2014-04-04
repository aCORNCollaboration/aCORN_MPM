#####################################################################
#
#  Name:         Makefile
#  Created by:   Michael Mendenhall
#
#  Contents:     Makefile for Michael's aCORN analysis code
#
#####################################################################


# assure correct shell is used
SHELL = /bin/sh
# apply implicit rules only for listed file types
.SUFFIXES:
.SUFFIXES: .c .cc .cpp .o
	 
# compiler command to use
CC = cc
CXX = g++

CXXFLAGS = -O3 -Wall -fPIC `root-config --cflags` -I. \
	-IIOUtils -IRootUtils -IBaseTypes -IMathUtils -ICalibration -IAnalysis -IStudies -IPhysics
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

VPATH = ./:IOUtils/:Physics/:ROOTUtils/

IOUtils = ControlMenu.o OutputManager.o PathUtils.o QFile.o SMExcept.o strutils.o

Physics = BetaSpectrum.o ElectronBindingEnergy.o FloatErr.o NuclEvtGen.o

ROOTUtils = GraphUtils.o

objects = $(IOUtils) $(Physics) $(ROOTUtils)


all: libaCORN_MPM.a

libaCORN_MPM.a: $(objects)
	ar rs libaCORN_MPM.a $(objects)

# generic rule for everything else .cc linked against libaCORN_MPM
% : %.cc libaCORN_MPM.a
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@


#aCORN Geant4:
#	mkdir -p g4build/
#	cd g4build; cmake -DGeant4_DIR=~/geant4.9.5/geant4.9.5-install/lib/Geant4-9.5.0/ ../ucnG4_dev/; make



#
# cleanup
#
.PHONY: clean
clean:
	-rm -f libaCORN_MPM.a
	-rm -f *.o
	-rm -rf *.dSYM

