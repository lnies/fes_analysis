# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

ARCH          = linux

# CXX           = 

CXX           = g++
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o 

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)


ifeq ($(ARCH),linux)
# Linux with egcs, gcc 2.9x, gcc 3.x (>= RedHat 5.2)
CXX           = g++
CXXFLAGS      = -Wall -g
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared
endif

# LDFlags was -O

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)-lMathMore -lSpectrum
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)-lMathMore -lSpectrum

#------------------------------------------------------------------------------

AW_PROTO60O   = aw_proto60.o aw_lib.o tagger_lib.o aio_module_def.o
AW_PROTO60S   = aw_proto60.$(SrcSuf) aw_lib.cpp tagger_lib.cpp aio_module_def.cpp
AW_PROTO60    = aw_proto60$(ExeSuf)

AW_MARKUSO   = aw_markus.o aw_lib.o tagger_lib.o aio_module_def.o
AW_MARKUSS   = aw_markus.$(SrcSuf) aw_lib.cpp tagger_lib.cpp aio_module_def.cpp
AW_MARKUS    = aw_markus$(ExeSuf)

AW_LUKASO   = aw_lukas.o aw_lib.o tagger_lib.o aio_module_def.o
AW_LUKASS   = aw_lukas.$(SrcSuf) aw_lib.cpp tagger_lib.cpp aio_module_def.cpp
AW_LUKAS    = aw_lukas$(ExeSuf)

AW_PROTO9O   = aw_proto9.o aw_lib.o tagger_lib.o aio_module_def.o
AW_PROTO9S   = aw_proto9.$(SrcSuf) aw_lib.cpp tagger_lib.cpp aio_module_def.cpp
AW_PROTO9    = aw_proto9$(ExeSuf)

AW_SHASHLIKO   = aw_shashlik.o aw_lib.o tagger_lib.o aio_module_def.o
AW_SHASHLIKS   = aw_shashlik.$(SrcSuf) aw_lib.cpp tagger_lib.cpp aio_module_def.cpp
AW_SHASHLIK    = aw_shashlik$(ExeSuf)

AW_TESTO   = aw_test.o aw_lib.o tagger_lib.o aio_module_def.o
AW_TESTS   = aw_test.$(SrcSuf) aw_lib.cpp tagger_lib.cpp aio_module_def.cpp
AW_TEST    = aw_test$(ExeSuf)

AW_ULIO   = aw_uli.o aw_lib_uli.o tagger_lib.o aio_module_def.o
AW_ULIS   = aw_uli.$(SrcSuf) aw_lib_uli.cpp tagger_lib.cpp aio_module_def.cpp
AW_ULI    = aw_uli$(ExeSuf)

OBJS          = $(AW_PROTO60O) 

PROGRAMS      = $(AW_PROTO60) 

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)


all:            $(PROGRAMS)

$(AW_PROTO60):     $(AW_PROTO60O)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"
$(AW_MARKUS):     $(AW_MARKUSO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"
$(AW_LUKAS):     $(AW_LUKASO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"
$(AW_PROTO9):     $(AW_PROTO9O)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"
$(AW_SHASHLIK):     $(AW_SHASHLIKO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"
$(AW_TEST):     $(AW_TESTO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"
$(AW_ULI):     $(AW_ULIO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"
aw_lib.o: aw_lib.cpp
		$(LD) -c aw_lib.cpp
aw_lib_uli.o: aw_lib_uli.cpp
		$(LD) -c aw_lib_uli.cpp

tagger_lib.o: tagger_lib.cpp
		$(LD) -c tagger_lib.cpp

aio_module_def.o: aio_module_def.cpp
		$(LD) -c aio_module_def.cpp

clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(PROGRAMS) $(EVENTSO) $(EVENTLIB) *Dict.* *.def *.exp \
		   *.root *.ps *.so .def so_locations
		@rm -rf cxx_repository

.SUFFIXES: .$(SrcSuf)

###


.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<

