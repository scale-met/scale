################################################################################
#
# Makefile for scale library
#
################################################################################

TOPDIR      = ../..
BUILD_DIR   = ./.libs
INCLUDE_DIR = ../include

SYSDEP_DIR = $(TOPDIR)/sysdep
DCUTILSDIR = $(TOPDIR)/dc_utils
GTOOLDIR   = $(TOPDIR)/gtool
LIBDIR     = $(TOPDIR)/lib
INCDIR     = $(TOPDIR)/include


LIBNAME = libscalelib.a

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)

VPATH = io:common:grid:atmos

VERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(VERSION),)
  VERSION = $(shell cat VERSION)
else
  VERSION := $(VERSION)
endif



OBJS =	\
	$(mod_rdma)	\
	mod_stdio.o	\
	mod_specfunc.o	\
	mod_process.o	\
	mod_const.o	\
	mod_calendar.o	\
	mod_random.o	\
	mod_papi.o	\
	mod_debug.o	\
	mod_time.o	\
	mod_misc.o	\
	mod_option.o	\
	mod_grid_cartesian.o	\
	mod_fileio.o	\
	mod_comm.o	\
	mod_topography.o	\
	mod_landuse.o	\
	mod_mapproj.o	\
	mod_grid_real.o	\
	mod_gridtrans.o	\
	mod_intepolation.o	\
	mod_stats.o	\
	mod_history.o	\
	mod_monitor.o	\
	mod_atmos_dyn_common.o	\
	mod_atmos_dyn_rk.o \
	mod_atmos_dyn_rk_heve.o \
	mod_atmos_dyn_rk_hevi.o \
	mod_atmos_dyn.o

all: makedir makedcutils makegtool $(LIBDIR)/$(LIBNAME) modules

makedcutils:
	cd $(DCUTILSDIR); $(MAKE)

makegtool:
	cd $(GTOOLDIR); $(MAKE)

makedir:
	mkdir -p $(BUILD_DIR)
	mkdir -p $(INCDIR)
	mkdir -p $(LIBDIR)

$(LIBDIR)/$(LIBNAME): $(BUILD_DIR)/$(LIBNAME)
	install $< $@

$(BUILD_DIR)/$(LIBNAME): $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

modules: 
	install $(patsubst %,$(BUILD_DIR)/%,$(MODS)) $(INCDIR)

allclean: distclean
	$(MAKE) -C $(DCUTILSDIR) allclean
	$(MAKE) -C $(GTOOLDIR) allclean
clean:
	rm -f $(PRJS) $(BUILD_DIR)/*.f90 $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.o $(BUILD_DIR)/*.lst

distclean: clean
	rm -f $(BUILD_DIR)/$(LIBNAME) $(BUILD_DIR)/*.a depend

depend:
	../tools/makedepend .

.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .c .erb .mod

%.F90 : %.F90.erb
	erb $< > $@

$(BUILD_DIR)/%.mod: %.F90
	make $(patsubst %.F90,$(BUILD_DIR)/%.o,$<)

$(INCDIR)/%.mod : $(BUILD_DIR)/%.mod
	install $< $@

$(BUILD_DIR)/%.o : %.F90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INCLUDE_DIR) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/mod_atmos_dyn.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INCLUDE_DIR) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/mod_atmos_dyn_rd.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INCDIR) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/mod_atmos_dyn_rd_heve.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INCDIR) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/mod_atmos_dyn_rd_hevi.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INCDIR) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/mod_atmos_phy_rd_profile.o:
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INCDIR) $(NETCDF_INCLUDE) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)

$(BUILD_DIR)/rdma.o					: rdma.c rdma.h

.PHONY : clean distclean allclean depend modules


-include depend
