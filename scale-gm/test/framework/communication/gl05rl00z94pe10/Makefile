################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..


# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
glevel      = 5
rlevel      = 0
nmpi        = 10
zlayer      = 94
vgrid       = vgrid94.dat

BINNAME     = communication

buildlocal:
	mkdir -p $(BUILD_DIR)
	$(MAKE) makelocal

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

BINS = $(BINDIR)/$(BINNAME)
LIBS = $(LIBDIR)/libscale.a $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a

ifneq ($(IGNOREDEP),T)
  DEPENDLIB = $(LIBS)
endif

VPATH = ..:                                       \
        $(BUILD_DIR):                             \
        $(SCALEGMDIR)/src/admin:                  \
        $(SCALEGMDIR)/src/atmos:                  \
        $(SCALEGMDIR)/src/atmos/dcmip:            \
        $(SCALEGMDIR)/src/atmos-gm/communication: \
        $(SCALEGMDIR)/src/atmos-gm/dynamics:      \
        $(SCALEGMDIR)/src/atmos-gm/grid:          \
        $(SCALEGMDIR)/src/atmos-physics/forcing:  \
        $(SCALEGMDIR)/src/atmos-physics/share:    \
        $(SCALEGMDIR)/src/common:                 \
        $(SCALEGMDIR)/src/io:                     \
        $(SCALEGMDIR)/src/preprocess:             \
        $(SCALEGMDIR)/src/tool:

VERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(VERSION),)
  VERSION  = $(shell cat VERSION)
else
  VERSION := $(VERSION)
endif

OBJS =	\
	mod_io_param.o \
	mod_adm.o      \
	mod_comm.o



makelocal:
	$(MAKE) makedir
	$(MAKE) -C $(DCUTILSDIR)
	$(MAKE) -C $(GTOOLDIR)
	$(MAKE) -C $(SCALELIBDIR)/src
	@echo;echo "Entering scale-gm...";echo "Current version is " $(VERSION)
	$(MAKE) makebin
	@echo "Complete making scale-gm."

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

makebin: $(BINS)

$(BINDIR)/$(BINNAME) : $(BUILD_DIR)/$(BINNAME).o $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(FC) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BUILD_DIR)/$(BINNAME).o : $(BINNAME).f90 $(patsubst %,$(BUILD_DIR)/%,$(OBJS))



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

%.mod: %.f90
	make $(patsubst %.f90,%.o,$<)

$(BUILD_DIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MODDIR) $(CONTRIB_INCLUDE) -I$(SCALELIBDIR)/include -I$(SCALEGMDIR)/include $(PAPI_INCLUDE) $(NETCDF_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)

$(BUILD_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -I$(SCALEGMDIR)/include -o $@ -c $<

.PHONY : clean distclean allclean dcutils gtool scalelib


ifneq ($(IGNOREDEP),T)
-include $(TOPDIR)/src/depend
endif
