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
index      = 500m40x40x40
precision  = real8
tracer     = kessler
sfcflx     = louis
turbulence = smg
radiation  = dummy
postfix    = 

VPATH = $(BUILD_DIR):include:common:io:grid:atmos:ocean:preprocess:user

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
	mod_stdio.o	\
	mod_process.o	\
	mod_const.o	\
	mod_calendar.o	\
	mod_random.o	\
	mod_specfunc.o	\
	mod_time.o	\
	mod_papi.o	\
	mod_misc.o	\
	mod_grid_cartesian_plane.o	\
	mod_fileio.o	\
	mod_geometrics.o	\
	$(mod_rdma)	\
	mod_comm.o	\
	mod_topography.o	\
	mod_history.o	\
	mod_monitor.o	\
	mod_atmos_sub_profile.o	\
	mod_atmos_sub_saturation.o	\
	mod_atmos_sub_thermodyn.o	\
	mod_atmos_sub_hydrostatic.o	\
	mod_atmos_vars.o	\
	mod_atmos_vars_sf.o	\
	mod_ocean_vars.o	\
	mod_atmos_sub_refstate.o	\
	mod_atmos_sub_boundary.o	\
	mod_atmos_sub_precipitation.o	\
	mod_atmos_dyn_fent_fullfct_lss.o	\
	mod_atmos_phy_sf_$(sfcflx).o	\
	mod_atmos_phy_tb_$(turbulence).o	\
	mod_atmos_phy_mp_common.o \
	mod_atmos_phy_mp_$(tracer).o \
	mod_atmos_phy_rd_$(radiation).o \
	mod_atmos.o	\
	mod_ocean_fixedsst.o	\
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

makebin: $(PRJS)
	install $(BUILD_DIR)/$(PRJ1) $(BINDIR)/$(BINNAME)
	install $(BUILD_DIR)/$(PRJ2) $(BINDIR)/$(INITNAME)

$(BUILD_DIR)/$(PRJ1): $(BUILD_DIR)/$(PRJ1).o $(LIBDIR)/$(LIBNAME) $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a
	$(FC) $(LFLAGS) -o $@ $^ $(NETCDF_LIBS) $(PAPI_LIBS)

$(BUILD_DIR)/$(PRJ2): $(BUILD_DIR)/$(PRJ2).o $(LIBDIR)/$(LIBNAME) $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a
	$(FC) $(LFLAGS) -o $@ $^ $(NETCDF_LIBS) $(PAPI_LIBS)

$(BUILD_DIR)/$(PRJ1).o: $(PRJ1).f90 $(LIBDIR)/$(LIBNAME)
$(BUILD_DIR)/$(PRJ2).o: $(PRJ2).f90 $(LIBDIR)/$(LIBNAME)

$(LIBDIR)/$(LIBNAME): $(BUILD_DIR)/$(LIBNAME)
	mkdir -p $(LIBDIR)
	install $< $@

$(BUILD_DIR)/$(LIBNAME): $(patsubst %,$(BUILD_DIR)/%,$(MODS))
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

$(BUILD_DIR)/mod_stdio.o	: mod_stdio.f90
$(BUILD_DIR)/mod_process.o	: mod_process.f90 mod_stdio.o
$(BUILD_DIR)/mod_const.o	: mod_const.f90 mod_stdio.o mod_process.o
$(BUILD_DIR)/mod_calendar.o	: mod_calendar.f90 mod_stdio.o mod_process.o
$(BUILD_DIR)/mod_random.o	: mod_random.f90 mod_stdio.o mod_process.o
$(BUILD_DIR)/mod_specfunc.o	: mod_specfunc.f90 mod_stdio.o mod_process.o
$(BUILD_DIR)/mod_time.o	: mod_time.f90 mod_stdio.o mod_process.o mod_calendar.o
$(BUILD_DIR)/mod_papi.o	: mod_papi.f90 mod_stdio.o mod_process.o
$(BUILD_DIR)/mod_misc.o	: mod_misc.f90 mod_stdio.o mod_process.o mod_time.o
$(BUILD_DIR)/mod_comm.o	: mod_comm.f90 mod_stdio.o mod_process.o mod_time.o mod_const.o mod_geometrics.o $(mod_rdma)
$(BUILD_DIR)/rdma.o	: rdma.c rdma.h

$(BUILD_DIR)/mod_grid_cartesian_plane.o	: mod_grid_cartesian_plane.f90 mod_stdio.o mod_process.o mod_time.o
$(BUILD_DIR)/mod_fileio.o	: mod_fileio.f90  mod_stdio.o mod_process.o mod_time.o mod_grid_cartesian_plane.o
$(BUILD_DIR)/mod_geometrics.o	: mod_geometrics.f90 mod_stdio.o mod_process.o mod_time.o mod_const.o mod_grid_cartesian_plane.o mod_fileio.o
$(BUILD_DIR)/mod_topography.o	: mod_topography.f90 mod_stdio.o mod_process.o mod_time.o mod_comm.o  mod_grid_cartesian_plane.o mod_fileio.o

$(BUILD_DIR)/mod_history.o	: mod_history.f90 mod_stdio.o mod_process.o mod_time.o mod_grid_cartesian_plane.o
$(BUILD_DIR)/mod_monitor.o	: mod_monitor.f90 mod_stdio.o mod_process.o mod_time.o mod_comm.o

$(BUILD_DIR)/mod_atmos.o	: mod_atmos.f90 mod_stdio.o mod_time.o mod_atmos_sub_thermodyn.o mod_atmos_sub_saturation.o mod_atmos_vars.o mod_atmos_sub_refstate.o mod_atmos_sub_boundary.o mod_atmos_dyn_fent_fullfct_lss.o mod_atmos_phy_sf_$(sfcflx).o mod_atmos_phy_tb_$(turbulence).o mod_atmos_phy_mp_$(tracer).o mod_atmos_phy_rd_$(radiation).o
$(BUILD_DIR)/mod_atmos_vars.o	: mod_atmos_vars.f90 mod_stdio.o mod_process.o mod_const.o mod_history.o mod_monitor.o mod_comm.o mod_grid_cartesian_plane.o mod_atmos_sub_thermodyn.o mod_time.o mod_misc.o mod_atmos_sub_saturation.o
$(BUILD_DIR)/mod_atmos_vars_sf.o	: mod_atmos_vars_sf.f90 mod_stdio.o mod_process.o mod_comm.o
$(BUILD_DIR)/mod_atmos_sub_profile.o	: mod_atmos_sub_profile.f90 mod_stdio.o mod_const.o mod_grid_cartesian_plane.o
$(BUILD_DIR)/mod_atmos_sub_thermodyn.o	: mod_atmos_sub_thermodyn.f90 mod_stdio.o mod_const.o
$(BUILD_DIR)/mod_atmos_sub_saturation.o	: mod_atmos_sub_saturation.f90 mod_stdio.o mod_const.o mod_process.o
$(BUILD_DIR)/mod_atmos_sub_hydrostatic.o	: mod_atmos_sub_hydrostatic.f90 mod_stdio.o mod_comm.o mod_process.o mod_const.o mod_grid_cartesian_plane.o
$(BUILD_DIR)/mod_atmos_sub_refstate.o	: mod_atmos_sub_refstate.f90 mod_stdio.o mod_process.o mod_const.o mod_grid_cartesian_plane.o mod_comm.o mod_fileio.o mod_atmos_sub_profile.o mod_atmos_sub_hydrostatic.o mod_atmos_sub_saturation.o mod_atmos_sub_thermodyn.o mod_atmos_vars.o
$(BUILD_DIR)/mod_atmos_sub_boundary.o	: mod_atmos_sub_boundary.f90 mod_stdio.o mod_process.o mod_const.o mod_grid_cartesian_plane.o mod_comm.o mod_time.o mod_atmos_vars.o
$(BUILD_DIR)/mod_atmos_sub_precipitation.o	: mod_atmos_sub_precipitation.f90 mod_stdio.o mod_time.o mod_const.o mod_grid_cartesian_plane.o mod_atmos_vars.o mod_atmos_sub_thermodyn.o

$(BUILD_DIR)/mod_atmos_dyn_fent_fullfct_lss.o	: mod_atmos_dyn_fent_fullfct_lss.f90 mod_stdio.o mod_time.o mod_process.o mod_grid_cartesian_plane.o mod_geometrics.o mod_atmos_vars.o mod_comm.o mod_const.o mod_history.o mod_atmos_sub_thermodyn.o mod_atmos_sub_refstate.o mod_atmos_sub_boundary.o mod_atmos_vars_sf.o mod_ocean_vars.o
$(BUILD_DIR)/mod_atmos_phy_sf_$(sfcflx).o	: mod_atmos_phy_sf_$(sfcflx).f90 mod_stdio.o mod_process.o mod_grid_cartesian_plane.o mod_const.o mod_atmos_vars.o mod_ocean_vars.o 
$(BUILD_DIR)/mod_atmos_phy_tb_$(turbulence).o	: mod_atmos_phy_tb_$(turbulence).f90 mod_stdio.o mod_grid_cartesian_plane.o mod_process.o mod_atmos_vars.o mod_comm.o mod_time.o mod_history.o mod_const.o
$(BUILD_DIR)/mod_atmos_phy_mp_common.o		: mod_atmos_phy_mp_common.f90 mod_stdio.o mod_const.o mod_process.o mod_time.o mod_grid_cartesian_plane.o mod_atmos_sub_thermodyn.o mod_atmos_sub_saturation.o
$(BUILD_DIR)/mod_atmos_phy_mp_$(tracer).o	: mod_atmos_phy_mp_$(tracer).f90 mod_stdio.o mod_const.o mod_time.o mod_process.o mod_grid_cartesian_plane.o mod_history.o mod_specfunc.o mod_atmos_vars.o mod_atmos_phy_mp_common.o mod_atmos_sub_precipitation.o mod_atmos_sub_thermodyn.o mod_atmos_sub_saturation.o
$(BUILD_DIR)/mod_atmos_phy_rd_$(radiation).o	: mod_atmos_phy_rd_$(radiation).f90 mod_stdio.o mod_process.o mod_comm.o mod_const.o mod_time.o mod_grid_cartesian_plane.o mod_history.o mod_atmos_vars.o 

$(BUILD_DIR)/mod_ocean.o	: mod_ocean.f90 mod_stdio.o mod_time.o mod_ocean_vars.o mod_ocean_fixedsst.o
$(BUILD_DIR)/mod_ocean_vars.o	: mod_ocean_vars.f90 mod_stdio.o mod_process.o mod_const.o mod_comm.o mod_time.o mod_grid_cartesian_plane.o
$(BUILD_DIR)/mod_ocean_fixedsst.o	: mod_ocean_fixedsst.f90 mod_stdio.o mod_process.o mod_ocean_vars.o mod_time.o mod_history.o

$(BUILD_DIR)/mod_mktopo.o	: mod_mktopo.f90 mod_stdio.o mod_process.o mod_grid_cartesian_plane.o mod_topography.o mod_const.o
$(BUILD_DIR)/mod_mkinit.o	: mod_mkinit.f90 mod_stdio.o mod_process.o mod_const.o mod_random.o mod_grid_cartesian_plane.o mod_atmos_vars.o mod_atmos_sub_profile.o mod_atmos_sub_hydrostatic.o mod_atmos_sub_saturation.o

$(BUILD_DIR)/mod_user.o	: mod_user.f90 mod_stdio.o mod_process.o mod_time.o mod_history.o mod_topography.o
$(BUILD_DIR)/scale3_init.o	: scale3_init.f90 mod_stdio.o mod_process.o mod_const.o mod_random.o mod_time.o mod_grid_cartesian_plane.o mod_geometrics.o mod_comm.o mod_topography.o mod_history.o mod_monitor.o mod_atmos_sub_thermodyn.o mod_atmos_sub_saturation.o mod_atmos_vars.o mod_mktopo.o mod_mkinit.o
$(BUILD_DIR)/scale.o	: scale.f90 mod_stdio.o mod_process.o mod_const.o mod_time.o mod_grid_cartesian_plane.o mod_geometrics.o mod_comm.o mod_topography.o mod_history.o mod_monitor.o mod_atmos.o mod_atmos_vars.o mod_ocean.o mod_ocean_vars.o

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
	$(FC) $(FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INCDIR) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/mod_atmos_dyn_fent_fullfct_lss.o:
	$(FC) $(FFLAGS_DYN) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INCDIR) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)


.PHONY : clean distclean allclean
