################################################################################
#
# Makefile for scale library
#
################################################################################

TOPDIR      = ../..
BUILD_DIR   = $(shell pwd)/.libs
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
        land:                       \
        land/grid:                  \
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

OBJS = $(scale_rdma)	\
	scale_atmos_dyn.o	\
	scale_atmos_dyn_common.o	\
	scale_atmos_dyn_rk.o	\
	scale_atmos_dyn_rk_heve.o	\
	scale_atmos_dyn_rk_hevi.o	\
	scale_atmos_dyn_rk_hivi.o	\
	scale_atmos_phy_ae.o	\
	scale_atmos_phy_ae_dummy.o	\
	scale_atmos_phy_mp.o	\
	scale_atmos_phy_mp_common.o	\
	scale_atmos_phy_mp_dry.o	\
	scale_atmos_phy_mp_kessler.o	\
	scale_atmos_phy_mp_sn14.o	\
	scale_atmos_phy_mp_suzuki10.o	\
	scale_atmos_phy_mp_sdm.o	\
	scale_atmos_phy_mp_tomita08.o	\
	scale_atmos_phy_rd.o	\
	scale_atmos_phy_rd_common.o	\
	scale_atmos_phy_rd_dummy.o	\
	scale_atmos_phy_rd_mstrnx.o	\
	scale_atmos_phy_rd_profile.o	\
	scale_atmos_phy_sf.o	\
	scale_atmos_phy_sf_const.o	\
	scale_atmos_phy_sf_louis.o	\
	scale_atmos_phy_tb.o	\
	scale_atmos_phy_tb_dummy.o	\
	scale_atmos_phy_tb_smg.o	\
	scale_atmos_phy_tb_dns.o	\
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
	scale_land_grid_cartesian.o	\
	scale_land_grid_index.o	\
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
	scale_tracer_sdm.o	\
	scale_tracer_suzuki10.o	\
	scale_tracer_tomita08.o

all:
	$(MAKE) makedir
	$(MAKE) makedcutils 
	$(MAKE) makegtool 
	$(MAKE) makesdm 
	$(MAKE) $(LIBDIR)/$(LIBNAME) 
	$(MAKE) modules

makedir:
	mkdir -p $(MODDIR)
	mkdir -p $(LIBDIR)
	mkdir -p $(BUILD_DIR)

makedcutils:
	$(MAKE) -C $(DCUTILSDIR)

makegtool:
	$(MAKE) -C $(GTOOLDIR)


ifeq ($(SDM),-D_SDM)
makesdm: $(BUILD_DIR)/scale_precision.o \
	$(BUILD_DIR)/scale_comm.o \
	$(BUILD_DIR)/scale_atmos_sub_thermodyn.o \
	$(BUILD_DIR)/scale_atmos_sub_saturation.o \
	$(BUILD_DIR)/scale_tracer.o \
	$(BUILD_DIR)/scale_stdio.o \
	$(BUILD_DIR)/scale_const.o \
	$(BUILD_DIR)/scale_grid_index.o \
	$(BUILD_DIR)/scale_topography.o \
	$(BUILD_DIR)/scale_grid_cartesian.o \
	$(BUILD_DIR)/scale_grid_real.o \
	$(BUILD_DIR)/scale_gridtrans.o \
	$(BUILD_DIR)/scale_process.o
	rsync -c $(SDMDIR)/scale_atmos_phy_mp_sdm.F90 $(SCALELIBDIR)/src/atmos-physics/microphysics/
	rsync -c $(SDMDIR)/scale_tracer_sdm.F90 $(SCALELIBDIR)/src/atmos-physics/microphysics/
	$(MAKE) -C $(MTDIR)
	$(MAKE) -C $(GAGETRYDIR)
	$(MAKE) -C $(SDMDIR)
else
makesdm:
	@echo "skip making additional component:SDM"
endif

$(LIBDIR)/$(LIBNAME): $(BUILD_DIR)/$(LIBNAME)
	install $< $@

$(BUILD_DIR)/$(LIBNAME): $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

modules:
	@install $(patsubst %,$(BUILD_DIR)/%,$(MODS)) $(MODDIR)



allclean: cleandcutils cleangtool cleangsdm distclean

cleandcutils:
	$(MAKE) -C $(DCUTILSDIR) allclean

cleangtool:
	$(MAKE) -C $(GTOOLDIR)   allclean

ifeq ($(SDM),-D_SDM)
cleangsdm:
	$(MAKE) -C $(GAGETRYDIR) allclean
	$(MAKE) -C $(MTDIR)      allclean
	$(MAKE) -C $(SDMDIR)     allclean
else
cleangsdm:
	@echo "skip cleaning additional component"
endif

distclean: clean
	rm -f $(BUILD_DIR)/$(LIBNAME) $(BUILD_DIR)/*.a

ifeq ($(SDM),-D_SDM)
clean:
	rm -f $(PRJS) $(BUILD_DIR)/*.f90 $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.o $(BUILD_DIR)/*.lst
	cp -pf $(SDMDIR)/scale_atmos_phy_mp_sdm_dummy.F90 $(SCALELIBDIR)/src/atmos-physics/microphysics/scale_atmos_phy_mp_sdm.F90
	cp -pf $(SDMDIR)/scale_tracer_sdm_dummy.F90 $(SCALELIBDIR)/src/atmos-physics/microphysics/scale_tracer_sdm.F90
else
clean:
	rm -f $(PRJS) $(BUILD_DIR)/*.f90 $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.o $(BUILD_DIR)/*.lst
endif

dependency:
	../tools/makedepend .



.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .c .erb .mod

%.F90 : %.F90.erb
	erb $< > $@

# Is this working? Could not build $(BUILD_DIR)/scale_precision.mod. \S.Shima.
$(BUILD_DIR)/%.mod: %.F90
	make $(patsubst %.F90,$(BUILD_DIR)/%.o,$<)

$(MODDIR)/%.mod : $(BUILD_DIR)/%.mod
	install $< $@

$(BUILD_DIR)/%.o : %.F90
	$(FC) $(FFLAGS)     $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) $(CONTRIB_INCLUDE) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/scale_atmos_dyn.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/scale_atmos_dyn_rk.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/scale_atmos_dyn_rk_heve.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/scale_atmos_dyn_rk_hevi.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/scale_atmos_dyn_rk_hivi.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/scale_atmos_phy_rd_profile.o:
	$(FC) $(FFLAGS)     $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) $(NETCDF_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
ifeq ($(SDM),-D_SDM)
$(BUILD_DIR)/scale_atmos_phy_mp_sdm.o: 
	$(FC) $(FFLAGS) $(FFLAGS_SDM_MONIT) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) $(CONTRIB_INCLUDE) -I$(SCALELIBDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
else
endif
$(BUILD_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -I$(SCALELIBDIR)/include -o $@ -c $<

$(BUILD_DIR)/scale_rdma.o : scale_rdma.c

.PHONY : clean distclean allclean depend modules dependency


-include depend
