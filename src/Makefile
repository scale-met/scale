################################################################################
#
# Makefile for main program
#
################################################################################

TOPDIR      = ..
BUILD_DIR   = ./.libs
INDEX_DIR   = ./include
SYSDEP_DIR  = ../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

# for build
index        = 500m40x40x40
precision    = real8
dynamics     = heve
tracer       = kessler
sfcflx       = louis
turbulence   = smg
microphysics = $(tracer)
radiation    = dummy
aerosol      = dummy
chemistry    = dummy
postfix      =

VPATH = $(BUILD_DIR):include:common:io:grid:atmos:land:ocean:cpl:preprocess:user

VERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(VERSION),)
  VERSION = $(shell cat VERSION)
else
  VERSION := $(VERSION)
endif



PRJ1 = scale3
PRJ2 = scale3_init

PRJS = $(BUILD_DIR)/$(PRJ1) $(BUILD_DIR)/$(PRJ2)

MODS =	\
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
	mod_grid_cartesian.o	\
	mod_fileio.o	\
	mod_comm.o	\
	mod_topography.o	\
	mod_landuse.o	\
	mod_grid_real.o	\
	mod_gridtrans.o	\
	mod_intepolation.o	\
	mod_stats.o	\
	mod_history.o	\
	mod_monitor.o	\
	mod_atmos_sub_profile.o	\
	mod_atmos_sub_saturation.o	\
	mod_atmos_sub_thermodyn.o	\
	mod_atmos_sub_hydrostatic.o	\
	mod_atmos_vars.o	\
	mod_atmos_vars_sf.o	\
	mod_land_vars.o	\
	mod_ocean_vars.o	\
	mod_cpl_vars.o \
	mod_cpl_atmos_land.o \
	mod_cpl.o \
	mod_atmos_sub_refstate.o	\
	mod_atmos_sub_boundary.o	\
	mod_atmos_sub_precipitation.o	\
	mod_atmos_sub_solarins.o	\
	mod_atmos_dyn_common.o	\
	mod_atmos_dyn_$(dynamics).o	\
	mod_atmos_dyn.o	\
	mod_atmos_phy_sf_$(sfcflx).o	\
	mod_atmos_phy_tb_$(turbulence).o	\
	mod_atmos_phy_mp_common.o \
	mod_atmos_phy_mp_$(microphysics).o \
	mod_atmos_phy_ae_$(aerosol).o \
	mod_atmos_phy_rd_common.o \
	mod_atmos_phy_rd_profile.o \
	mod_atmos_phy_rd_$(radiation).o \
	mod_atmos.o	\
	mod_land_phy_bucket.o	\
	mod_land.o	\
	mod_ocean_phy_fixedsst.o	\
	mod_ocean.o	\
	mod_mktopo.o	\
	mod_mkinit.o	\
	mod_user.o

all: makedir makedcutils makegtool makeinclude makebin

makedcutils:
	cd $(DCUTILSDIR); $(MAKE)

makegtool:
	cd $(GTOOLDIR); $(MAKE)

makedir:
	mkdir -p $(BINDIR)
	mkdir -p $(BUILD_DIR)

makeinclude:
	cp $(INDEX_DIR)/inc_index_$(index).f90      $(BUILD_DIR)/inc_index.h
	cp $(INCDIR)/inc_tracer_$(tracer).f90       $(BUILD_DIR)/inc_tracer.h
	cp $(INCDIR)/inc_precision_$(precision).f90 $(BUILD_DIR)/inc_precision.h
	cp $(INCDIR)/inc_land.f90                   $(BUILD_DIR)/inc_land.h

makebin: $(PRJS)
	install $(BUILD_DIR)/$(PRJ1) $(BINDIR)/$(BINNAME)
	install $(BUILD_DIR)/$(PRJ2) $(BINDIR)/$(INITNAME)

$(BUILD_DIR)/$(PRJ1): $(BUILD_DIR)/$(PRJ1).o $(LIBDIR)/$(LIBNAME) $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a
	$(FC) $(LDFLAGS) -o $@ $^ $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BUILD_DIR)/$(PRJ2): $(BUILD_DIR)/$(PRJ2).o $(LIBDIR)/$(LIBNAME) $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a
	$(FC) $(LDFLAGS) -o $@ $^ $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BUILD_DIR)/$(PRJ1).o: $(PRJ1).f90 $(LIBDIR)/$(LIBNAME)
$(BUILD_DIR)/$(PRJ2).o: $(PRJ2).f90 $(LIBDIR)/$(LIBNAME)

$(LIBDIR)/$(LIBNAME): $(BUILD_DIR)/$(LIBNAME)
	mkdir -p $(LIBDIR)
	install $< $@

$(BUILD_DIR)/$(LIBNAME): $(patsubst %,$(BUILD_DIR)/%,$(MODS))
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

$(BUILD_DIR)/rdma.o					: rdma.c rdma.h

$(BUILD_DIR)/mod_specfunc.o			: mod_specfunc.f90
$(BUILD_DIR)/mod_process.o				: mod_process.f90				mod_stdio.o
$(BUILD_DIR)/mod_const.o				: mod_const.f90				mod_stdio.o mod_process.o
$(BUILD_DIR)/mod_calendar.o			: mod_calendar.f90			mod_stdio.o mod_process.o
$(BUILD_DIR)/mod_random.o				: mod_random.f90				mod_stdio.o mod_process.o
$(BUILD_DIR)/mod_papi.o					: mod_papi.f90					mod_stdio.o mod_process.o
$(BUILD_DIR)/mod_debug.o				: mod_debug.f90				mod_const.o
$(BUILD_DIR)/mod_time.o					: mod_time.f90					mod_stdio.o mod_process.o mod_const.o mod_calendar.o
$(BUILD_DIR)/mod_misc.o					: mod_misc.f90					mod_stdio.o mod_process.o mod_time.o
$(BUILD_DIR)/mod_comm.o					: mod_comm.f90					mod_stdio.o mod_process.o mod_time.o $(mod_rdma)
$(BUILD_DIR)/mod_stats.o				: mod_stats.f90				mod_stdio.o mod_process.o mod_time.o mod_comm.o mod_grid_real.o

$(BUILD_DIR)/mod_grid_cartesian.o	: mod_grid_cartesian.f90	mod_stdio.o mod_process.o mod_time.o
$(BUILD_DIR)/mod_topography.o			: mod_topography.f90			mod_stdio.o mod_process.o mod_time.o mod_comm.o mod_fileio.o
$(BUILD_DIR)/mod_landuse.o				: mod_landuse.f90			   mod_stdio.o mod_process.o mod_time.o mod_comm.o mod_fileio.o
$(BUILD_DIR)/mod_grid_real.o			: mod_grid_real.f90			mod_stdio.o mod_process.o mod_time.o mod_const.o mod_fileio.o mod_grid_cartesian.o mod_topography.o
$(BUILD_DIR)/mod_gridtrans.o			: mod_gridtrans.f90			mod_stdio.o mod_process.o mod_time.o mod_comm.o mod_grid_cartesian.o mod_grid_real.o
$(BUILD_DIR)/mod_intepolation.o		: mod_intepolation.f90		mod_stdio.o mod_process.o mod_time.o mod_const.o mod_grid_cartesian.o mod_grid_real.o

$(BUILD_DIR)/mod_stdio.o				: mod_stdio.f90
$(BUILD_DIR)/mod_fileio.o				: mod_fileio.f90				mod_stdio.o mod_process.o mod_time.o mod_grid_cartesian.o
$(BUILD_DIR)/mod_history.o				: mod_history.f90		   	mod_stdio.o mod_process.o mod_time.o mod_grid_cartesian.o mod_intepolation.o
$(BUILD_DIR)/mod_monitor.o				: mod_monitor.f90			   mod_stdio.o mod_process.o mod_time.o mod_stats.o

$(BUILD_DIR)/mod_atmos.o	: mod_atmos.f90 mod_stdio.o mod_time.o mod_atmos_sub_thermodyn.o mod_atmos_sub_saturation.o mod_atmos_vars.o mod_atmos_sub_refstate.o mod_atmos_sub_boundary.o \
                             mod_atmos_dyn.o mod_atmos_phy_sf_$(sfcflx).o mod_atmos_phy_tb_$(turbulence).o mod_atmos_phy_mp_$(microphysics).o mod_atmos_phy_rd_$(radiation).o mod_cpl_vars.o
$(BUILD_DIR)/mod_atmos_vars.o	: mod_atmos_vars.f90 mod_stdio.o mod_process.o mod_const.o mod_comm.o mod_stats.o mod_history.o mod_monitor.o mod_grid_cartesian.o mod_atmos_sub_thermodyn.o mod_time.o mod_misc.o mod_atmos_sub_saturation.o
$(BUILD_DIR)/mod_atmos_vars_sf.o	: mod_atmos_vars_sf.f90 mod_stdio.o mod_process.o mod_comm.o
$(BUILD_DIR)/mod_atmos_sub_profile.o	: mod_atmos_sub_profile.f90 mod_stdio.o mod_const.o mod_grid_cartesian.o
$(BUILD_DIR)/mod_atmos_sub_thermodyn.o	: mod_atmos_sub_thermodyn.f90 mod_stdio.o mod_const.o
$(BUILD_DIR)/mod_atmos_sub_saturation.o	: mod_atmos_sub_saturation.f90 mod_stdio.o mod_const.o mod_process.o
$(BUILD_DIR)/mod_atmos_sub_hydrostatic.o	: mod_atmos_sub_hydrostatic.f90 mod_stdio.o mod_comm.o mod_process.o mod_const.o mod_grid_cartesian.o
$(BUILD_DIR)/mod_atmos_sub_refstate.o	: mod_atmos_sub_refstate.f90 mod_stdio.o mod_process.o mod_const.o mod_grid_cartesian.o mod_comm.o mod_fileio.o mod_atmos_sub_profile.o mod_atmos_sub_hydrostatic.o mod_atmos_sub_saturation.o mod_atmos_sub_thermodyn.o mod_atmos_vars.o
$(BUILD_DIR)/mod_atmos_sub_boundary.o	: mod_atmos_sub_boundary.f90 mod_stdio.o mod_process.o mod_const.o mod_grid_cartesian.o mod_comm.o mod_time.o mod_atmos_vars.o
$(BUILD_DIR)/mod_atmos_sub_precipitation.o	: mod_atmos_sub_precipitation.f90 mod_stdio.o mod_time.o mod_const.o mod_grid_cartesian.o mod_atmos_vars.o mod_atmos_sub_thermodyn.o
$(BUILD_DIR)/mod_atmos_sub_solarins.o	: mod_atmos_sub_solarins.f90 mod_stdio.o mod_const.o mod_calendar.o

$(BUILD_DIR)/mod_atmos_dyn_common.o	: mod_atmos_dyn_common.f90 mod_comm.o mod_debug.o
$(BUILD_DIR)/mod_atmos_dyn_$(dynamics).o	: mod_atmos_dyn_$(dynamics).f90 mod_comm.o mod_const.o mod_debug.o mod_atmos_dyn_common.o mod_atmos_vars.o mod_process.o
$(BUILD_DIR)/mod_atmos_dyn.o	: mod_atmos_dyn.f90 mod_atmos_dyn_$(dynamics).o mod_atmos_dyn_common.o mod_stdio.o mod_time.o mod_process.o mod_grid_cartesian.o mod_grid_real.o mod_gridtrans.o mod_atmos_vars.o mod_comm.o mod_const.o mod_history.o mod_atmos_sub_thermodyn.o mod_atmos_sub_refstate.o mod_atmos_sub_boundary.o mod_atmos_vars_sf.o mod_ocean_vars.o mod_debug.o
$(BUILD_DIR)/mod_atmos_phy_sf_$(sfcflx).o	: mod_atmos_phy_sf_$(sfcflx).f90 mod_stdio.o mod_process.o mod_grid_cartesian.o mod_const.o mod_atmos_vars.o mod_atmos_vars_sf.o mod_ocean_vars.o
$(BUILD_DIR)/mod_atmos_phy_tb_$(turbulence).o	: mod_atmos_phy_tb_$(turbulence).f90 mod_stdio.o mod_grid_cartesian.o mod_process.o mod_atmos_vars.o mod_comm.o mod_time.o mod_history.o mod_const.o mod_debug.o
$(BUILD_DIR)/mod_atmos_phy_mp_common.o		: mod_atmos_phy_mp_common.f90 mod_stdio.o mod_const.o mod_process.o mod_time.o mod_grid_cartesian.o mod_atmos_sub_thermodyn.o mod_atmos_sub_saturation.o
$(BUILD_DIR)/mod_atmos_phy_mp_$(microphysics).o	: mod_atmos_phy_mp_$(microphysics).f90 mod_stdio.o mod_const.o mod_time.o mod_process.o mod_specfunc.o mod_grid_cartesian.o mod_history.o mod_atmos_vars.o mod_atmos_phy_mp_common.o mod_atmos_sub_precipitation.o mod_atmos_sub_thermodyn.o mod_atmos_sub_saturation.o
$(BUILD_DIR)/mod_atmos_phy_ae_$(aerosol).o	: mod_atmos_phy_ae_$(aerosol).f90
$(BUILD_DIR)/mod_atmos_phy_rd_common.o		: mod_atmos_phy_rd_common.f90 mod_stdio.o mod_time.o mod_grid_cartesian.o
$(BUILD_DIR)/mod_atmos_phy_rd_$(radiation).o	: mod_atmos_phy_rd_$(radiation).f90 mod_stdio.o mod_process.o mod_comm.o mod_const.o mod_time.o mod_grid_cartesian.o mod_history.o mod_atmos_sub_solarins.o mod_atmos_phy_rd_profile.o mod_atmos_phy_mp_$(tracer).o mod_atmos_phy_ae_$(aerosol).o mod_atmos_phy_rd_common.o mod_atmos_vars.o
$(BUILD_DIR)/mod_atmos_phy_rd_profile.o	: mod_atmos_phy_rd_profile.f90 mod_stdio.o mod_process.o mod_const.o mod_grid_cartesian.o mod_calendar.o

$(BUILD_DIR)/mod_land.o	: mod_land.f90 mod_stdio.o mod_time.o mod_land_vars.o mod_land_phy_bucket.o mod_cpl_vars.o
$(BUILD_DIR)/mod_land_vars.o	: mod_land_vars.f90 mod_stdio.o mod_process.o mod_const.o mod_comm.o mod_time.o mod_misc.o mod_fileio.o
$(BUILD_DIR)/mod_land_phy_bucket.o	: mod_land_phy_bucket.f90 mod_stdio.o mod_process.o mod_time.o mod_land_vars.o

$(BUILD_DIR)/mod_ocean.o	: mod_ocean.f90 mod_stdio.o mod_time.o mod_ocean_vars.o mod_ocean_phy_fixedsst.o
$(BUILD_DIR)/mod_ocean_vars.o	: mod_ocean_vars.f90 mod_stdio.o mod_process.o mod_const.o mod_comm.o mod_time.o mod_misc.o mod_fileio.o
$(BUILD_DIR)/mod_ocean_phy_fixedsst.o	: mod_ocean_phy_fixedsst.f90 mod_stdio.o mod_process.o mod_time.o mod_ocean_vars.o

$(BUILD_DIR)/mod_cpl.o	: mod_cpl.f90 mod_stdio.o mod_time.o mod_cpl_vars.o mod_cpl_atmos_land.o
$(BUILD_DIR)/mod_cpl_vars.o	: mod_cpl_vars.f90 mod_stdio.o mod_process.o mod_const.o mod_comm.o mod_time.o mod_misc.o mod_fileio.o
$(BUILD_DIR)/mod_cpl_atmos_land.o	: mod_cpl_atmos_land.f90 mod_stdio.o mod_const.o mod_grid_cartesian.o mod_cpl_vars.o

$(BUILD_DIR)/mod_mktopo.o	: mod_mktopo.f90 mod_stdio.o mod_process.o mod_grid_cartesian.o mod_topography.o mod_const.o
$(BUILD_DIR)/mod_mkinit.o	: mod_mkinit.f90 mod_stdio.o mod_process.o mod_const.o mod_random.o mod_grid_cartesian.o mod_atmos_vars.o mod_atmos_sub_profile.o mod_atmos_sub_hydrostatic.o mod_atmos_sub_saturation.o

$(BUILD_DIR)/mod_user.o	: mod_user.f90 mod_stdio.o mod_process.o mod_time.o mod_history.o

$(BUILD_DIR)/scale3_init.o	: scale3_init.f90	mod_stdio.o mod_process.o mod_const.o mod_random.o mod_time.o mod_fileio.o mod_grid_cartesian.o mod_comm.o mod_topography.o mod_landuse.o mod_grid_real.o mod_gridtrans.o mod_intepolation.o mod_stats.o mod_history.o mod_atmos_sub_thermodyn.o mod_atmos_sub_saturation.o mod_atmos_vars.o mod_mktopo.o mod_mkinit.o
$(BUILD_DIR)/scale.o			: scale.f90			mod_stdio.o mod_process.o mod_const.o mod_random.o mod_time.o mod_fileio.o mod_grid_cartesian.o mod_comm.o mod_topography.o mod_landuse.o mod_grid_real.o mod_gridtrans.o mod_intepolation.o mod_stats.o mod_history.o mod_monitor.o mod_atmos.o mod_atmos_vars.o mod_ocean.o mod_ocean_vars.o mod_cpl.o mod_cpl_vars.o

allclean: distclean
	$(MAKE) -C $(DCUTILSDIR) allclean
	$(MAKE) -C $(GTOOLDIR) allclean
clean:
	rm -f $(PRJS) $(BUILD_DIR)/*.h $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.o $(BUILD_DIR)/*.lst

distclean: clean
	rm -f $(BUILD_DIR)/$(LIBNAME) $(BUILD_DIR)/*.a

.SUFFIXES:
.SUFFIXES: .o .f90 .c

$(BUILD_DIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INCDIR) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/mod_atmos_dyn.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INCDIR) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/mod_atmos_dyn_$(dynamics).o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INCDIR) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/mod_atmos_phy_rd_profile.o:
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INCDIR) $(NETCDF_INCLUDE) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)


.PHONY : clean distclean allclean
