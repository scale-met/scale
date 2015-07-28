################################################################################
#
# Makefile for scale library
#
################################################################################

TOPDIR      = ../..
SYSDEP_DIR := ../../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

LIBNAME = libscale.a

DEPENDLIB = $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a

VPATH = common:                     \
        common/tracer:              \
        coupler:                    \
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
        atmos-les/io:               \
        atmos-les/dynamics:         \
        atmos-les/forcing:          \
        land:                       \
        land/grid:                  \
        ocean:                      \
        urban:                      \
        urban/grid:

LIBVERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(LIBVERSION),)
  LIBVERSION  = $(shell cat VERSION)
else
  LIBVERSION := $(LIBVERSION)
endif

ifeq ($(RDMA),-D_USE_RDMA)
  scale_rdma := scale_rdma.o
else
  scale_rdma :=
endif

ifeq ($(MP),)
  OBJS_MP = scale_atmos_phy_mp_dry.o	\
            scale_atmos_phy_mp_kessler.o	\
            scale_atmos_phy_mp_sn14.o	\
            scale_atmos_phy_mp_suzuki10.o	\
            scale_atmos_phy_mp_tomita08.o
else
  OBJS_MP = scale_atmos_phy_mp_$(MP).o
  FFLAGS += -DMP=$(MP)
endif

ifeq ($(AE),)
  OBJS_AE = scale_atmos_phy_ae_dummy.o \
	    scale_atmos_phy_ae_kajino13.o
else
  OBJS_AE = scale_atmos_phy_ae_$(AE).o
endif

ifeq ($(RD),)
  OBJS_RD = scale_atmos_phy_rd_mstrnx.o	\
            scale_atmos_phy_rd_mm5sw.o	\
            scale_atmos_phy_rd_profile.o
else
  OBJS_RD = scale_atmos_phy_rd_$(RD).o
  FFLAGS += -DRD=$(RD)
endif

ifeq ($(TB),)
  OBJS_TB = scale_atmos_phy_tb_dummy.o	\
            scale_atmos_phy_tb_smg.o	\
            scale_atmos_phy_tb_dns.o	\
            scale_atmos_phy_tb_mynn.o	\
            scale_atmos_phy_tb_hybrid.o
else
  OBJS_TB = scale_atmos_phy_tb_$(TB).o
  FFLAGS += -DTB=$(TB)
endif

ifeq ($(TRACER),)
  OBJS_TRACER = scale_tracer_dry.o	\
                scale_tracer_kessler.o	\
                scale_tracer_sn14.o	\
                scale_tracer_suzuki10.o	\
                scale_tracer_tomita08.o
else
  OBJS_TRACER = scale_tracer_$(TRACER).o
  FFLAGS += -DTRACER=$(TRACER)
endif

ifeq ($(AETRACER),)
  OBJS_AETRACER = scale_aetracer_none.o	\
                scale_aetracer_kajino13.o
else
  OBJS_AETRACER = scale_aetracer_$(AETRACER).o
endif

OBJS_NAME = $(scale_rdma)	\
	scale_atmos_dyn.o	\
	scale_atmos_dyn_common.o	\
	scale_atmos_dyn_rk.o	\
	scale_atmos_dyn_rk_heve.o	\
	scale_atmos_dyn_rk_hevi.o	\
	scale_atmos_dyn_rk_hivi.o	\
	scale_atmos_phy_ae.o	\
	scale_atmos_phy_ae_dummy.o	\
	$(OBJS_AE)	\
	scale_atmos_phy_mp.o	\
	scale_atmos_phy_mp_common.o	\
	$(OBJS_MP)	\
	scale_atmos_phy_rd.o	\
	scale_atmos_phy_rd_common.o	\
	$(OBJS_RD)	\
	scale_atmos_phy_sf.o	\
	scale_atmos_phy_sf_bulkcoef.o	\
	scale_atmos_phy_sf_bulk.o	\
	scale_atmos_phy_sf_const.o	\
	scale_atmos_phy_tb.o	\
	scale_atmos_phy_tb_common.o	\
	$(OBJS_TB)	\
	scale_atmos_sub_bottom.o	\
	scale_atmos_sub_boundary.o	\
	scale_atmos_sub_hydrostatic.o	\
	scale_atmos_sub_profile.o	\
	scale_atmos_sub_refstate.o	\
	scale_atmos_sub_saturation.o	\
	scale_atmos_sub_solarins.o	\
	scale_atmos_sub_thermodyn.o	\
	scale_bulkflux.o	\
	scale_calendar.o	\
	scale_comm.o	\
	scale_const.o	\
	scale_debug.o	\
	scale_external_io.o	\
	scale_fileio.o	\
	scale_grid_cartesian.o	\
	scale_grid_index.o	\
	scale_grid_nest.o	\
	scale_grid_real.o	\
	scale_gridtrans.o	\
	scale_history.o	\
	scale_intepolation.o	\
	scale_interpolation_nest.o	\
	scale_land_grid.o	\
	scale_land_grid_index.o	\
	scale_land_phy.o	\
	scale_land_phy_slab.o	\
	scale_land_phy_matsiro.o	\
	scale_land_sfc.o	\
	scale_land_sfc_slab.o	\
	scale_landuse.o	\
	scale_les_process.o	\
	scale_mapproj.o	\
	scale_matrix.o	\
	scale_monitor.o	\
	scale_ocean_phy.o	\
	scale_ocean_phy_slab.o	\
	scale_ocean_sfc.o	\
	scale_ocean_sfc_slab.o	\
	scale_precision.o	\
	scale_process.o	\
	scale_prof.o	\
	scale_random.o	\
	scale_roughness.o	\
	scale_specfunc.o	\
	scale_statistics.o	\
	scale_stdio.o	\
	scale_time.o	\
	scale_topography.o	\
	scale_tracer.o	\
	$(OBJS_TRACER)	\
	$(OBJS_AETRACER)	\
	scale_urban_grid.o	\
	scale_urban_grid_index.o	\
	scale_urban_phy.o	\
	scale_urban_phy_slc.o

OBJS =  $(patsubst %,$(OBJ_DIR)/%,$(OBJS_NAME))


all:
	$(MAKE) makedir
	@echo;echo "Entering scalelib..."; echo "Current version is " $(LIBVERSION)
	$(MAKE) -C $(DCUTILSDIR)
	$(MAKE) -C $(GTOOLDIR)
	$(MAKE) $(LIBDIR)/$(LIBNAME)
	$(MAKE) modules
	@echo "Complete making scalelib."

makedir:
	@mkdir -p $(MODDIR)
	mkdir -p $(LIBDIR)

$(LIBDIR)/libdcutils.a:
	$(MAKE) -C $(DCUTILSDIR)

$(LIBDIR)/libgtool.a: $(LIBDIR)/libdcutils.a
	$(MAKE) -C $(GTOOLDIR)

$(LIBNAME): $(LIBDIR)/$(LIBNAME)
	install $< $@

$(LIBDIR)/$(LIBNAME): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

ifeq ($(IGNOREDEP),)
modules: $(OBJS)
	install $(MODS) $(MODDIR)
else
modules:
endif


allclean: distclean
	$(MAKE) -C $(DCUTILSDIR) allclean
	$(MAKE) -C $(GTOOLDIR)   allclean
	rm -f $(LIBDIR)/$(LIBNAME)

distclean: clean
	rm -f $(LIBNAME)

clean:
	rm -f *.o *.mod *.lst

dependency:
	./makedepend .



.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .c .erb .mod

%.F90 : %.F90.erb
	erb $< > $@

%.mod: %.F90
	make $(patsubst %.F90,%.o,$<)

$(MODDIR)/%.mod : %.mod
	install $< $@

$(OBJ_DIR)/%.o : %.F90
	$(FC) $(FFLAGS)     $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(LIBVERSION)\" -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(INDEX_DIR) -I$(SCALELIBDIR)/include -I$(OBJ_DIR) $(PAPI_INCLUDE) $(NETCDF_INCLUDE) $(MODDIROPT) $(OBJ_DIR) -o $@ -c $<
$(OBJ_DIR)/scale_atmos_dyn.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(LIBVERSION)\" -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(INDEX_DIR) -I$(SCALELIBDIR)/include -I$(OBJ_DIR) $(PAPI_INCLUDE) $(MODDIROPT) $(OBJ_DIR) -o $@ -c $<
$(OBJ_DIR)/scale_atmos_dyn_common.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(LIBVERSION)\" -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(INDEX_DIR) -I$(SCALELIBDIR)/include -I$(OBJ_DIR) $(PAPI_INCLUDE) $(MODDIROPT) $(OBJ_DIR) -o $@ -c $<
$(OBJ_DIR)/scale_atmos_dyn_rk.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(LIBVERSION)\" -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(INDEX_DIR) -I$(SCALELIBDIR)/include -I$(OBJ_DIR) $(PAPI_INCLUDE) $(MODDIROPT) $(OBJ_DIR) -o $@ -c $<
$(OBJ_DIR)/scale_atmos_dyn_rk_heve.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(LIBVERSION)\" -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(INDEX_DIR) -I$(SCALELIBDIR)/include -I$(OBJ_DIR) $(PAPI_INCLUDE) $(MODDIROPT) $(OBJ_DIR) -o $@ -c $<
$(OBJ_DIR)/scale_atmos_dyn_rk_hevi.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(LIBVERSION)\" -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(INDEX_DIR) -I$(SCALELIBDIR)/include -I$(OBJ_DIR) $(PAPI_INCLUDE) $(MODDIROPT) $(OBJ_DIR) -o $@ -c $<
$(OBJ_DIR)/scale_atmos_dyn_rk_hivi.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(LIBVERSION)\" -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(INDEX_DIR) -I$(SCALELIBDIR)/include -I$(OBJ_DIR) $(PAPI_INCLUDE) $(MODDIROPT) $(OBJ_DIR) -o $@ -c $<

$(OBJ_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -I$(SCALELIBDIR)/include -o $@ -c $<

$(OBJ_DIR)/scale_rdma.o : scale_rdma.c

.PHONY : clean distclean allclean depend modules dependency


ifeq ($(IGNOREDEP),)
-include depend
endif
