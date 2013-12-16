################################################################################
#
# Makefile for scale library
#
################################################################################

TOPDIR      = ../..
BUILD_DIR   = ./.libs
INCLUDE_DIR = ../include

include $(TOPDIR)/Mkinclude
include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)

LIBNAME = libscalelib.a

VPATH = io:common:grid:atmos

VERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(VERSION),)
  VERSION = $(shell cat VERSION)
else
  VERSION := $(VERSION)
endif

ifeq ($(RDMA),-D_USE_RDMA)
  mod_rdma := rdma.o
else
  mod_rdma :=
endif


OBJS =	\
	$(mod_rdma)	\
	mod_precision.o	\
	mod_index.o	\
	mod_tracer.o	\
	mod_tracer_dry.o	\
	mod_tracer_kessler.o	\
	mod_tracer_tomita08.o	\
	mod_tracer_sn13.o	\
	mod_tracer_sn13w.o	\
	mod_tracer_binw.o	\
	mod_tracer_binf.o	\
	mod_tracer_hbinw.o	\
	mod_tracer_hbinf.o	\
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
	mod_atmos_sub_saturation.o	\
	mod_atmos_sub_thermodyn.o	\
	mod_atmos_sub_solarins.o	\
	mod_atmos_dyn_common.o	\
	mod_atmos_dyn_rk.o \
	mod_atmos_dyn_rk_heve.o \
	mod_atmos_dyn_rk_hevi.o \
	mod_atmos_dyn.o	\
	mod_atmos_phy_sf_const.o	\
	mod_atmos_phy_sf_louis.o	\
	mod_atmos_phy_sf.o	\
	mod_atmos_phy_tb_smg.o	\
	mod_atmos_phy_tb_dummy.o	\
	mod_atmos_phy_tb.o	\
	mod_atmos_phy_mp_common.o	\
	mod_atmos_phy_mp_dry.o	\
	mod_atmos_phy_mp_kessler.o	\
	mod_atmos_phy_mp_tomita08.o	\
	mod_atmos_phy_mp_sn13.o	\
	mod_atmos_phy_mp_sn13w.o	\
	mod_atmos_phy_mp_binw.o	\
	mod_atmos_phy_mp_binf.o	\
	mod_atmos_phy_mp_hbinw.o	\
	mod_atmos_phy_mp_hbinf.o	\
	mod_atmos_phy_mp.o	\
	mod_atmos_phy_rd_common.o	\
	mod_atmos_phy_rd_profile.o	\
	mod_atmos_phy_rd_dummy.o	\
	mod_atmos_phy_rd_dycoms2.o	\
	mod_atmos_phy_rd_mstrnX.o	\
	mod_atmos_phy_rd.o	\
	mod_atmos_phy_ae_dummy.o	\
	mod_atmos_phy_ae.o

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
	@install $(patsubst %,$(BUILD_DIR)/%,$(MODS)) $(INCDIR)

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
