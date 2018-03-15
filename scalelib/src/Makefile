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

DEPENDLIB =

VPATH = common:                          \
        coupler:                         \
        io:                              \
	file:                            \
	grid:                            \
	communication:                   \
        atmosphere/grid:                 \
        atmosphere/physics/pbl:          \
        atmosphere/physics/microphysics: \
        atmosphere/physics/cumulus:      \
        atmosphere/physics/radiation:    \
        atmosphere/physics/surface:      \
        atmosphere/physics/aerosol:      \
        atmosphere/physics/chemistry:    \
        atmos-physics:                   \
        atmos-physics/share:             \
        atmos-physics/microphysics:      \
        atmos-physics/surfaceflux:       \
        atmos-physics/turbulence:        \
        atmos-rm:                        \
        atmos-rm/communication:          \
        atmos-rm/grid:                   \
        atmos-rm/io:                     \
        atmos-rm/numeric:                \
        atmos-rm/dynamics:               \
        atmos-rm/forcing:                \
        land:                            \
        land/grid:                       \
        ocean:                           \
        ocean/grid:                      \
        urban:                           \
        urban/grid:

LIBVERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(LIBVERSION),)
  LIBVERSION  = $(shell cat VERSION)
else
  LIBVERSION := $(LIBVERSION)
endif

SDM_EXIST := $(shell find $(TOPDIR)/contrib -name SDM 2> /dev/null )

OBJS_DYN = scale_atmos_dyn_fvm_flux_cd2.o	\
           scale_atmos_dyn_fvm_flux_cd4.o	\
           scale_atmos_dyn_fvm_flux_cd6.o	\
           scale_atmos_dyn_fvm_flux_ud1.o	\
           scale_atmos_dyn_fvm_flux_ud3.o	\
           scale_atmos_dyn_fvm_flux_ud3Koren1993.o	\
           scale_atmos_dyn_fvm_flux_ud5.o	\
           scale_atmos_dyn_fvm_flux.o	\
           scale_atmos_dyn_tinteg_large_euler.o	\
           scale_atmos_dyn_tinteg_large_rk3.o	\
           scale_atmos_dyn_tinteg_large.o	\
           scale_atmos_dyn_tinteg_short_rk3.o	\
           scale_atmos_dyn_tinteg_short_rk4.o	\
           scale_atmos_dyn_tinteg_short.o	\
           scale_atmos_dyn_tinteg_tracer_euler.o	\
           scale_atmos_dyn_tinteg_tracer_rk3.o	\
           scale_atmos_dyn_tinteg_tracer.o	\
           scale_atmos_dyn_tstep_large_fvm_heve.o	\
           scale_atmos_dyn_tstep_large.o	\
           scale_atmos_dyn_tstep_short_fvm_heve.o	\
           scale_atmos_dyn_tstep_short_fvm_hevi.o	\
           scale_atmos_dyn_tstep_short_fvm_hivi.o	\
           scale_atmos_dyn_tstep_short.o	\
           scale_atmos_dyn_tstep_tracer_fvm_heve.o	\
           scale_atmos_dyn_tstep_tracer.o

OBJS_MP = scale_atmos_phy_mp_dry.o	\
          scale_atmos_phy_mp_kessler.o	\
          scale_atmos_phy_mp_sn14.o	\
          scale_atmos_phy_mp_sdm.o	\
          scale_atmos_phy_mp_suzuki10.o	\
          scale_atmos_phy_mp_tomita08.o

OBJS_AE = scale_atmos_phy_ae_kajino13.o

OBJS_CH = scale_atmos_phy_ch_rn222.o

OBJS_RD = scale_atmos_phy_rd_mstrnx.o	\
          scale_atmos_phy_rd_offline.o	\
          scale_atmos_phy_rd_profile.o

OBJS_SF = scale_atmos_phy_sf_bulk.o	\
          scale_atmos_phy_sf_const.o

OBJS_TB = scale_atmos_phy_tb_smg.o	\
          scale_atmos_phy_tb_d1980.o	\
          scale_atmos_phy_tb_dns.o

OBJS_BL = scale_atmos_phy_bl_mynn.o

OBJS_CP = scale_atmos_phy_cp_kf.o



OBJS_NAME = \
	scale.o \
	scale_file_netcdf.o	\
	scale_file_f.o	\
	scale_file_h.o	\
	scale_file.o	\
	scale_file_cartesC.o	\
	scale_file_history.o	\
	scale_file_history_cartesC.o	\
	scale_file_external_input.o	\
	scale_file_external_input_cartesC.o	\
	scale_atmos_grid_cartesC_index.o	\
	scale_atmos_grid_cartesC.o	\
	scale_atmos_grid_cartesC_real.o	\
	scale_atmos_grid_cartesC_metric.o	\
	scale_atmos_dyn.o	\
	scale_atmos_dyn_common.o	\
	$(OBJS_DYN)	\
	$(OBJS_NUMERIC_FDM)	\
	\
	scale_atmos_phy_mp.o	\
	scale_atmos_phy_mp_common.o	\
	scale_atmos_phy_mp_convert.o	\
	$(OBJS_MP)	\
	\
	$(OBJS_AE)	\
	\
	$(OBJS_CH)	\
	\
	scale_atmos_phy_rd_common.o	\
	$(OBJS_RD)	\
	\
	$(OBJS_SF)	\
	\
	scale_atmos_phy_tb.o	\
	scale_atmos_phy_tb_common.o	\
	$(OBJS_TB)	\
        \
	$(OBJS_BL)	\
	\
	scale_atmos_phy_cp_common.o	\
	$(OBJS_CP)	\
	\
	scale_atmos_sub_bottom.o	\
	scale_atmos_sub_boundary.o	\
	scale_atmos_sub_hydrostatic.o	\
	scale_atmos_sub_profile.o	\
	scale_atmos_sub_refstate.o	\
	scale_atmos_sub_saturation.o	\
	scale_atmos_sub_solarins.o	\
	scale_atmos_sub_thermodyn.o	\
	scale_atmos_sub_aerosol.o	\
	scale_atmos_sub_hydrometeor.o	\
	scale_atmos_sub_adiabat.o	\
	scale_atmos_sub_diagnostic.o	\
	scale_bulkflux.o	\
	scale_calendar.o	\
	scale_comm_cartesC.o	\
	scale_comm_cartesC_nest.o	\
	scale_const.o	\
	scale_debug.o	\
	scale_index.o		\
	scale_interp.o	\
	scale_interp_vert.o	\
	scale_ocean_grid_cartesC_index.o	\
	scale_ocean_grid_cartesC.o	\
	scale_ocean_grid_cartesC_real.o	\
	scale_ocean_phy.o	\
	scale_ocean_phy_const.o	\
	scale_ocean_phy_slab.o	\
	scale_ocean_sfc.o	\
	scale_ocean_sfc_const.o	\
	scale_ocean_sfc_slab.o	\
	scale_ocean_sfc.o	\
	scale_ocean_sfc_const.o	\
	scale_ocean_sfc_slab.o	\
	scale_land_grid_cartesC_index.o	\
	scale_land_grid_cartesC.o	\
	scale_land_grid_cartesC_real.o	\
	scale_land_phy.o	\
	scale_land_phy_const.o	\
	scale_land_phy_slab.o	\
	scale_land_phy_matsiro.o	\
	scale_land_phy_snow_ky90.o	\
	scale_land_phy_snow_diagnos.o	\
	scale_land_sfc.o	\
	scale_land_sfc_const.o	\
	scale_land_sfc_thin-slab.o	\
	scale_land_sfc_thick-slab.o	\
	scale_urban_grid_cartesC_index.o	\
	scale_urban_grid_cartesC.o	\
	scale_urban_grid_cartesC_real.o	\
	scale_urban_phy.o	\
	scale_urban_phy_slc.o	\
	scale_landuse.o	\
	scale_rm_process.o	\
	scale_fpm.o  \
	scale_mapprojection.o	\
	scale_matrix.o	\
	scale_statistics.o	\
	scale_monitor.o	\
	scale_monitor_cartesC.o	\
	scale_precision.o	\
	scale_process.o	\
	scale_prof.o	\
	scale_random.o	\
	scale_roughness.o	\
	scale_sigvars.o \
	scale_specfunc.o	\
	scale_stdio.o	\
	scale_sys_sigvars.o  \
	scale_time.o	\
	scale_topography.o	\
	scale_tracer.o	\
	scale_vector.o

OBJS = $(patsubst %,$(OBJ_DIR)/%,$(OBJS_NAME))

OBJS_FOR_DYN = $(OBJ_DIR)/scale_atmos_dyn.o	\
	$(OBJ_DIR)/scale_atmos_dyn_common.o	\
	$(patsubst %,$(OBJ_DIR)/%,$(OBJS_DYN))

OBJS_FOR_SDM = $(OBJ_DIR)/scale_precision.o	\
	$(OBJ_DIR)/scale_prof.o	\
	$(OBJ_DIR)/scale_comm.o	\
	$(OBJ_DIR)/scale_atmos_sub_thermodyn.o	\
	$(OBJ_DIR)/scale_atmos_sub_saturation.o	\
	$(OBJ_DIR)/scale_atmos_phy_mp_common.o	\
	$(OBJ_DIR)/scale_tracer.o	\
	$(OBJ_DIR)/scale_tracer_sdm.o	\
	$(OBJ_DIR)/scale_stdio.o	\
	$(OBJ_DIR)/scale_history.o	\
	$(OBJ_DIR)/scale_const.o	\
	$(OBJ_DIR)/scale_time.o	\
	$(OBJ_DIR)/scale_atmos_grid_cartesC_index.o	\
	$(OBJ_DIR)/scale_topography.o	\
	$(OBJ_DIR)/scale_grid_cartesian.o	\
	$(OBJ_DIR)/scale_atmos_grid_cartesC_real.o	\
	$(OBJ_DIR)/scale_gridtrans.o	\
	$(OBJ_DIR)/scale_rm_process.o

all:
	$(MAKE) makedir
	@echo;echo "Entering scalelib..."; echo "Current version is " $(LIBVERSION)
ifeq ($(SCALE_ENABLE_SDM),T)
	$(MAKE) checksdm
	$(MAKE) overridesdm
	$(MAKE) makesdm
endif
	$(MAKE) $(LIBDIR)/$(LIBNAME)
	$(MAKE) modules
	@echo "Complete making scalelib."

makedir:
	mkdir -p $(OBJ_DIR)
	mkdir -p $(MODDIR)
	mkdir -p $(LIBDIR)

$(LIBNAME): $(LIBDIR)/$(LIBNAME)
	install $< $@

$(LIBDIR)/$(LIBNAME): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

ifeq ($(SCALE_IGNORE_SRCDEP),T)
modules:
else
modules: $(OBJS)
	install $(MODS) $(MODDIR)
endif

makesdm: $(OBJS_FOR_SDM)
	$(MAKE) -C $(MTDIR)
	$(MAKE) -C $(GAGETRYDIR)
	$(MAKE) -C $(SDMDIR)

checksdm:
ifneq ($(SDM_EXIST),$(SDMDIR))
	@echo "SDM directory is missing!" $(SDM_EXIST)
	@echo "********************************************************************************************* "
	@echo "** SDM module, which is compiled with 'SCALE_ENABLE_SDM=T', is not supported.              ** "
	@echo "** If you want to use SDM module, please contact Dr. S.Shima(s_shima@sim.u-hyogo.ac.jp)    ** "
	@echo "********************************************************************************************* "
	exit 1
endif

overridesdm:
	@echo;echo "Copy SDM-related file from $(SDMDIR)"
	rsync -th $(SDMDIR)/scale_atmos_phy_mp_sdm.F90 $(SCALELIBDIR)/src/atmos-physics/microphysics/
	rsync -th $(SDMDIR)/scale_tracer_sdm.F90       $(SCALELIBDIR)/src/common/tracer/
	@echo

restoresdm:
	@echo;echo "Restore SDM-related file"
	rsync -th $(SDMDIR)/scale_atmos_phy_mp_sdm_dummy.F90 $(SCALELIBDIR)/src/atmos-physics/microphysics/scale_atmos_phy_mp_sdm.F90
	rsync -th $(SDMDIR)/scale_tracer_sdm_dummy.F90       $(SCALELIBDIR)/src/common/tracer/scale_tracer_sdm.F90
	@echo

allcleansdm:
	$(MAKE) -C $(MTDIR)      allclean
	$(MAKE) -C $(GAGETRYDIR) allclean
	$(MAKE) -C $(SDMDIR)     allclean

cleansdm:
	$(MAKE) -C $(SDMDIR)     clean



allclean: distclean
ifeq ($(SCALE_ENABLE_SDM),T)
	$(MAKE) checksdm
	$(MAKE) restoresdm
	$(MAKE) allcleansdm
endif
	rm -f $(LIBDIR)/$(LIBNAME)

distclean: clean
	rm -f $(LIBNAME)

clean:
ifeq ($(SCALE_ENABLE_SDM),T)
	$(MAKE) checksdm
	$(MAKE) restoresdm
	$(MAKE) cleansdm
endif
	rm -rf .libs
	rm -f *.o *.mod *.lst

dependency:
	./makedepend .



.SUFFIXES:
.SUFFIXES: .o .F90 .c .erb .mod

%.F90 : %.F90.erb
	erb $< > $@

%.mod: %.F90
	make $(patsubst %.F90,%.o,$<)

$(MODDIR)/%.mod : %.mod
	install $< $@

$(OBJ_DIR)/%.o : %.F90
	$(FC) $(FFLAGS)     $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(LIBVERSION)\" -I$(SCALELIBDIR)/include -I$(OBJ_DIR) $(SCALE_PAPI_INCLUDE) $(SCALE_NETCDF_INCLUDE) $(CONTRIB_INCLUDE) $(MODDIROPT) $(OBJ_DIR) -o $@ -c $<
$(OBJS_FOR_DYN) :
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(LIBVERSION)\" -I$(SCALELIBDIR)/include -I$(OBJ_DIR) $(SCALE_PAPI_INCLUDE) $(MODDIROPT) $(OBJ_DIR) -o $@ -c $<

file/scale_file.F90 : scale_file.F90.erb

file/scale_file_history.F90 : scale_file_history.F90.erb

scale_atmos_dyn_fvm_flux_%.F90 : scale_atmos_dyn_fvm_flux_udcd.F90.erb
	fname=$@ erb $< > $@
#	erb fname=$@ $< > $@

$(OBJ_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -I$(SCALELIBDIR)/include $(SCALE_NETCDF_INCLUDE) -o $@ -c $<

.PHONY : clean distclean allclean depend modules dependency

include $(TOPDIR)/utils/make/Make_environments

ifneq ($(SCALE_IGNORE_SRCDEP),T)
-include depend
endif
# DO NOT DELETE
