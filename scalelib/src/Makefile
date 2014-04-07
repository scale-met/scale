################################################################################
#
# Makefile for scale library
#
################################################################################

TOPDIR      = ../..
BUILD_DIR   = ./.libs
SYSDEP_DIR := ../../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

LIBNAME = libscale.a

VPATH = common:                     \
        io:                         \
        atmos-physics:              \
        atmos-physics/share:        \
        atmos-physics/cumulus:      \
        atmos-physics/microphysics: \
        atmos-physics/radiation:    \
        atmos-physics/surfaceflux:  \
        atmos-physics/turbulence:   \
        atmos-physics/aerosol:      \
        atmos-physics/chemistry:    \
        atmos-les:                  \
        atmos-les/communication:    \
        atmos-les/grid:             \
        atmos-les/dynamics:         \
        ocean:                      \
        cpl:

VERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(VERSION),)
  VERSION = $(shell cat VERSION)
else
  VERSION := $(VERSION)
endif

ifeq ($(RDMA),-D_USE_RDMA)
  scale_rdma := scale_rdma.o
else
  scale_rdma :=
endif

OBJS = $(mod_rdma)	\
	scale_atmos_dyn.o	\
	scale_atmos_dyn_common.o	\
	scale_atmos_dyn_rk.o	\
	scale_atmos_dyn_rk_heve.o	\
	scale_atmos_dyn_rk_hevi.o	\
	scale_atmos_phy_ae.o	\
	scale_atmos_phy_ae_dummy.o	\
	scale_atmos_phy_mp.o	\
	scale_atmos_phy_mp_common.o	\
	scale_atmos_phy_mp_dry.o	\
	scale_atmos_phy_mp_kessler.o	\
	scale_atmos_phy_mp_sn14.o	\
	scale_atmos_phy_mp_suzuki10.o	\
	scale_atmos_phy_mp_tomita08.o	\
	scale_atmos_phy_rd.o	\
	scale_atmos_phy_rd_common.o	\
	scale_atmos_phy_rd_dummy.o	\
	scale_atmos_phy_rd_dycoms2.o	\
	scale_atmos_phy_rd_mstrnx.o	\
	scale_atmos_phy_rd_profile.o	\
	scale_atmos_phy_sf.o	\
	scale_atmos_phy_sf_const.o	\
	scale_atmos_phy_sf_louis.o	\
	scale_atmos_phy_sf_rico.o	\
	scale_atmos_phy_tb.o	\
	scale_atmos_phy_tb_dummy.o	\
	scale_atmos_phy_tb_smg.o	\
	scale_atmos_sub_boundary.o	\
	scale_atmos_sub_hydrostatic.o	\
	scale_atmos_sub_profile.o	\
	scale_atmos_sub_refstate.o	\
	scale_atmos_sub_saturation.o	\
	scale_atmos_sub_solarins.o	\
	scale_atmos_sub_thermodyn.o	\
	scale_calendar.o	\
	scale_comm.o	\
	scale_const.o	\
	scale_cpl_atmos_land.o	\
	scale_cpl_atmos_land_bulk.o	\
	scale_cpl_atmos_land_const.o	\
	scale_cpl_atmos_ocean.o	\
	scale_cpl_atmos_ocean_bulk.o	\
	scale_cpl_atmos_ocean_const.o	\
	scale_cpl_sub_bulkcoef.o	\
	scale_debug.o	\
	scale_fileio.o	\
	scale_grid_cartesian.o	\
	scale_grid_index.o	\
	scale_grid_real.o	\
	scale_gridtrans.o	\
	scale_history.o	\
	scale_intepolation.o	\
	scale_landuse.o	\
	scale_mapproj.o	\
	scale_monitor.o	\
	scale_ocean_sub_roughness.o	\
	scale_precision.o	\
	scale_process.o	\
	scale_prof.o	\
	scale_random.o	\
	scale_specfunc.o	\
	scale_stats.o	\
	scale_stdio.o	\
	scale_time.o	\
	scale_topography.o	\
	scale_tracer.o	\
	scale_tracer_dry.o	\
	scale_tracer_kessler.o	\
	scale_tracer_sn14.o	\
	scale_tracer_suzuki10.o	\
	scale_tracer_tomita08.o

all: makedir makedcutils makegtool $(LIBDIR)/$(LIBNAME) modules

makedir:
	mkdir -p $(BUILD_DIR)
	mkdir -p $(MODDIR)
	mkdir -p $(LIBDIR)

makedcutils:
	$(MAKE) -C $(DCUTILSDIR)

makegtool:
	$(MAKE) -C $(GTOOLDIR)

$(LIBDIR)/$(LIBNAME): $(BUILD_DIR)/$(LIBNAME)
	install $< $@

$(BUILD_DIR)/$(LIBNAME): $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

modules:
	@install $(patsubst %,$(BUILD_DIR)/%,$(MODS)) $(MODDIR)



allclean: cleandcutils cleangtool distclean

cleandcutils:
	$(MAKE) -C $(DCUTILSDIR) allclean

cleangtool:
	$(MAKE) -C $(GTOOLDIR)   allclean

distclean: clean
	rm -f $(BUILD_DIR)/$(LIBNAME) $(BUILD_DIR)/*.a depend

clean:
	rm -f $(PRJS) $(BUILD_DIR)/*.f90 $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.o $(BUILD_DIR)/*.lst

depend:
	../tools/makedepend .



.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .c .erb .mod

%.F90 : %.F90.erb
	erb $< > $@

$(BUILD_DIR)/%.mod: %.F90
	make $(patsubst %.F90,$(BUILD_DIR)/%.o,$<)

$(MODDIR)/%.mod : $(BUILD_DIR)/%.mod
	install $< $@

$(BUILD_DIR)/%.o : %.F90
	$(FC) $(FFLAGS)     $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/scale_atmos_dyn.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/scale_atmos_dyn_common.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/scale_atmos_dyn_rk.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/scale_atmos_dyn_rk_heve.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/scale_atmos_dyn_rk_hevi.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/scale_atmos_phy_rd_profile.o:
	$(FC) $(FFLAGS)     $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) $(NETCDF_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)

$(BUILD_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)

$(BUILD_DIR)/scale_rdma.o : scale_rdma.c scale_rdma.h

.PHONY : clean distclean allclean depend modules


-include depend
