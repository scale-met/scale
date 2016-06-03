################################################################################
#
# Makefile for main program
#
################################################################################

TOPDIR      = $(abspath ../../../../..)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = $(TOPDIR)/sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

PPNAME   = NONE
INITNAME = NONE
BINNAME  = scale-rm_fileio

DEPENDLIB = $(LIBDIR)/libscale.a $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a

VERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(VERSION),)
  VERSION  = $(shell cat VERSION)
else
  VERSION := $(VERSION)
endif

OBJS =	\
	mod_ioapi.o



all:
	$(MAKE) makedir
	$(MAKE) -C $(DCUTILSDIR)
	$(MAKE) -C $(GTOOLDIR)
	$(MAKE) -C $(SCALELIBDIR)/src
	@echo;echo "Entering scale-rm...";echo "Current version is " $(VERSION)
	$(MAKE) makebin
	@echo "Complete making scale-rm."

makedir:
	mkdir -p $(BINDIR)
	mkdir -p $(LIBDIR)
	mkdir -p $(BUILD_DIR)

$(LIBDIR)/libdcutils.a :
	$(MAKE) -C $(DCUTILSDIR)

$(LIBDIR)/libgtool.a : $(LIBDIR)/libdcutils.a
	$(MAKE) -C $(GTOOLDIR)

$(LIBDIR)/libscale.a : $(LIBDIR)/libdcutils.a $(LIBDIR)/libgtool.a
	$(MAKE) -C $(SCALELIBDIR)/src

makebin: $(BINDIR)/$(BINNAME)

$(BINDIR)/$(BINNAME)  : $(BUILD_DIR)/$(BINNAME).o  $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o
	$(FC) $(LDFLAGS) -o $@ $^ $(DEPENDLIB) $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BUILD_DIR)/$(BINNAME).o  : $(BINNAME).f90  $(DEPENDLIB) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o



allclean: distclean
	$(MAKE) -C $(DCUTILSDIR)      allclean
	$(MAKE) -C $(GTOOLDIR)        allclean
	$(MAKE) -C $(SCALELIBDIR)/src allclean
	rm -f $(BINDIR)/$(BINNAME)

distclean: clean
	rm -f $(BUILD_DIR)/$(BINNAME)

clean:
	rm -f $(BUILD_DIR)/*.*



.SUFFIXES:
.SUFFIXES: .o .f90 .c .mod

$(BUILD_DIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MODDIR) -I$(SCALELIBDIR)/include -I$(SCALERMDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)

%.mod: %.f90
	make $(patsubst %.f90,%.o,$<)


.PHONY: clean dcutils gtool scalelib

# parameters for run
PPCONF      = NONE
INITCONF    = NONE
RUNCONF     = sample.conf
TPROC       = 6

# required data (parameters,distributed files)
DATDIR      = ../../../data
DATPARAM    =
DATDISTS    =

run: jobshell
	$(JOBSUB) run.sh

jobshell:
	@if [ -f ../Mkjobshell.$(SCALE_SYS).sh ]; then \
		sh ../Mkjobshell.$(SCALE_SYS).sh $(BINDIR) $(INITNAME) $(BINNAME) $(INITCONF) $(RUNCONF) $(TPROC) \
		                                 $(DATDIR) $(DATPARAM) $(DATDISTS); \
	else \
		sh $(SYSDEP_DIR)/Mkjobshell.$(SCALE_SYS).sh $(BINDIR) $(INITNAME) $(BINNAME) $(INITCONF) $(RUNCONF) $(TPROC) \
		                                            $(DATDIR) $(DATPARAM) $(DATDISTS); \
	fi


