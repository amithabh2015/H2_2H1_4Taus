ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++

#CXXFLAGS      = -g -ggdb3 -O2 -pipe -Wall -Wextra -fPIC -ansi -D_GNU_SOURCE  
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
CXXFLAGS     += $(ROOTCFLAGS)

NGLIBS        = $(ROOTGLIBS) 
#NGLIBS        += -lMinuit -lFoam -lASImage
#NGLIBS        += -lMinuit -lRooFit -lRooFitCore -lFoam -lASImage

LD            = g++
LDFLAGS       = -g
OFLAGS        = -fPIC -c $(ROOTCFLAGS)
SOFLAGS       = -fPIC -shared $(ROOTCFLAGS)

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)

GLIBS         = $(NGLIBS)

INCLUDEDIR    = ./interface/
SRCDIR        = ./src/
CXX	     += -I$(INCLUDEDIR)

OUTLIB	      = $(CMSSW_BASE)/lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: $(CMSSW_BASE)/lib/


all:  DecayAnglesCalculator

DecayAnglesCalculator: $(SRCDIR)DecayAnglesCalculator.cc $(INCLUDEDIR)DecayAnglesCalculator.h
# 	$(CXX) $(CXXFLAGS) -o DecayAnglesCalculator $(GLIBS) $ DecayAnglesCalculator.cc
	$(CXX) $(OFLAGS) $(GLIBS) $(SRCDIR)DecayAnglesCalculator.cc
	$(LD) $(SOFLAGS) $(GLIBS) -o DecayAnglesCalculator.so DecayAnglesCalculator.o

clean:
	rm -f $(OUTLIB)*.o $(OUTLIBCOMMON)*.o
	rm -f $(OUTLIB)*.so
