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
        io:                              \
        file:                            \
        communication:                   \
        atmosphere/common:               \
        atmosphere/grid:                 \
        atmosphere/physics/aerosol:      \
        atmosphere/physics/chemistry:    \
        atmosphere/physics/cumulus:      \
        atmosphere/physics/microphysics: \
        atmosphere/physics/pbl:          \
        atmosphere/physics/radiation:    \
        atmosphere/physics/surface:      \
        land:                            \
        land/grid:                       \
        land/dynamics:                   \
        land/physics:                    \
        ocean:                           \
        ocean/grid:                      \
        ocean/dynamics:                  \
        ocean/physics:                   \
        urban:                           \
        urban/grid:                      \
        urban/dynamics:                  \
        coupler:                         \
        coupler/physics/surface:         \
        \
        atmos-rm/dynamics:               \
        atmos-rm/turbulence:

LIBVERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(LIBVERSION),)
  LIBVERSION  = $(shell cat VERSION)
else
  LIBVERSION := $(LIBVERSION)
endif

SDM_EXIST := $(shell find $(TOPDIR)/contrib -name SDM 2> /dev/null )

OBJS = \
	$(OBJ_DIR)/scale.o	\
	\
	$(OBJ_DIR)/scale_sys_sigvars.o	\
	$(OBJ_DIR)/scale_file_netcdf.o	\
	$(OBJ_DIR)/scale_file_f.o	\
	\
	$(OBJ_DIR)/scale_calendar.o	\
	$(OBJ_DIR)/scale_const.o	\
	$(OBJ_DIR)/scale_debug.o	\
	$(OBJ_DIR)/scale_fpm.o	\
	$(OBJ_DIR)/scale_index.o	\
	$(OBJ_DIR)/scale_interp.o	\
	$(OBJ_DIR)/scale_interp_vert.o	\
	$(OBJ_DIR)/scale_mapprojection.o	\
	$(OBJ_DIR)/scale_coriolis.o	\
	$(OBJ_DIR)/scale_matrix.o	\
	$(OBJ_DIR)/scale_sort.o	\
	$(OBJ_DIR)/scale_filter.o	\
	$(OBJ_DIR)/scale_prc.o	\
	$(OBJ_DIR)/scale_prc_cartesC.o	\
	$(OBJ_DIR)/scale_prc_icoA.o	\
	$(OBJ_DIR)/scale_precision.o	\
	$(OBJ_DIR)/scale_prof.o	\
	$(OBJ_DIR)/scale_random.o	\
	$(OBJ_DIR)/scale_sigvars.o	\
	$(OBJ_DIR)/scale_specfunc.o	\
	$(OBJ_DIR)/scale_statistics.o	\
	$(OBJ_DIR)/scale_time.o	\
	$(OBJ_DIR)/scale_topography.o	\
	$(OBJ_DIR)/scale_tracer.o	\
	$(OBJ_DIR)/scale_vector.o	\
	$(OBJ_DIR)/scale_io.o	\
	$(OBJ_DIR)/scale_monitor.o	\
	$(OBJ_DIR)/scale_monitor_cartesC.o	\
	$(OBJ_DIR)/scale_comm_cartesC.o	\
	$(OBJ_DIR)/scale_comm_cartesC_nest.o	\
	$(OBJ_DIR)/scale_comm_icoA.o	\
	$(OBJ_DIR)/scale_file.o	\
	$(OBJ_DIR)/scale_file_cartesC.o	\
	$(OBJ_DIR)/scale_file_external_input.o	\
	$(OBJ_DIR)/scale_file_external_input_cartesC.o	\
	$(OBJ_DIR)/scale_file_h.o	\
	$(OBJ_DIR)/scale_file_history.o	\
	$(OBJ_DIR)/scale_file_history_cartesC.o	\
	$(OBJ_DIR)/scale_file_tiledata.o	\
	\
	$(OBJ_DIR)/scale_atmos_adiabat.o	\
	$(OBJ_DIR)/scale_atmos_aerosol.o	\
	$(OBJ_DIR)/scale_atmos_bottom.o	\
	$(OBJ_DIR)/scale_atmos_diagnostic.o	\
	$(OBJ_DIR)/scale_atmos_diagnostic_cartesC.o	\
	$(OBJ_DIR)/scale_atmos_hydrometeor.o	\
	$(OBJ_DIR)/scale_atmos_hydrostatic.o	\
	$(OBJ_DIR)/scale_atmos_profile.o	\
	$(OBJ_DIR)/scale_atmos_refstate.o	\
	$(OBJ_DIR)/scale_atmos_saturation.o	\
	$(OBJ_DIR)/scale_atmos_solarins.o	\
	$(OBJ_DIR)/scale_atmos_thermodyn.o	\
	$(OBJ_DIR)/scale_atmos_grid_cartesC.o	\
	$(OBJ_DIR)/scale_atmos_grid_cartesC_index.o	\
	$(OBJ_DIR)/scale_atmos_grid_cartesC_metric.o	\
	$(OBJ_DIR)/scale_atmos_grid_cartesC_real.o	\
	$(OBJ_DIR)/scale_atmos_grid_icoA.o	\
	$(OBJ_DIR)/scale_atmos_grid_icoA_index.o	\
	$(OBJ_DIR)/scale_atmos_dyn.o	\
	$(OBJ_DIR)/scale_atmos_dyn_common.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux_cd2.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux_cd4.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux_cd6.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux_ud1.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux_ud3.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux_ud3Koren1993.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux_ud5.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_large.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_large_euler.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_large_rk3.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_short.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_short_rk3.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_short_rk4.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_tracer.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_tracer_euler.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_tracer_rk3.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_large.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_large_fvm_heve.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_short.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_short_fvm_heve.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_short_fvm_hevi.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_short_fvm_hivi.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_tracer.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_tracer_fvm_heve.o	\
	$(OBJ_DIR)/scale_atmos_phy_ae_kajino13.o	\
	$(OBJ_DIR)/scale_atmos_phy_bl_mynn.o	\
	$(OBJ_DIR)/scale_atmos_phy_ch_rn222.o	\
	$(OBJ_DIR)/scale_atmos_sfc_ch_rn222.o	\
	$(OBJ_DIR)/scale_atmos_phy_cp_common.o	\
	$(OBJ_DIR)/scale_atmos_phy_cp_kf.o	\
	$(OBJ_DIR)/scale_atmos_phy_mp_common.o	\
	$(OBJ_DIR)/scale_atmos_phy_mp_kessler.o	\
	$(OBJ_DIR)/scale_atmos_phy_mp_sn14.o	\
	$(OBJ_DIR)/scale_atmos_phy_mp_suzuki10.o	\
	$(OBJ_DIR)/scale_atmos_phy_mp_tomita08.o	\
	$(OBJ_DIR)/scale_atmos_phy_rd_common.o	\
	$(OBJ_DIR)/scale_atmos_phy_rd_mstrnx.o	\
	$(OBJ_DIR)/scale_atmos_phy_rd_offline.o	\
	$(OBJ_DIR)/scale_atmos_phy_rd_profile.o	\
	$(OBJ_DIR)/scale_atmos_phy_sf_bulk.o	\
	$(OBJ_DIR)/scale_atmos_phy_sf_const.o	\
	\
	$(OBJ_DIR)/scale_atmos_phy_tb_common.o	\
	$(OBJ_DIR)/scale_atmos_phy_tb_d1980.o	\
	$(OBJ_DIR)/scale_atmos_phy_tb_dns.o	\
	$(OBJ_DIR)/scale_atmos_phy_tb_smg.o	\
	\
	$(OBJ_DIR)/scale_ocean_grid_cartesC.o	\
	$(OBJ_DIR)/scale_ocean_grid_cartesC_index.o	\
	$(OBJ_DIR)/scale_ocean_grid_cartesC_real.o	\
	$(OBJ_DIR)/scale_ocean_grid_icoA.o	\
	$(OBJ_DIR)/scale_ocean_grid_icoA_index.o	\
	$(OBJ_DIR)/scale_ocean_dyn_slab.o	\
	$(OBJ_DIR)/scale_ocean_dyn_offline.o	\
	$(OBJ_DIR)/scale_ocean_phy_ice_simple.o	\
	$(OBJ_DIR)/scale_ocean_phy_albedo.o	\
	$(OBJ_DIR)/scale_ocean_phy_albedo_nakajima00.o	\
	$(OBJ_DIR)/scale_ocean_phy_roughness.o	\
	$(OBJ_DIR)/scale_ocean_phy_roughness_miller92.o	\
	$(OBJ_DIR)/scale_ocean_phy_roughness_moon07.o	\
	$(OBJ_DIR)/scale_ocean_phy_tc.o	\
	$(OBJ_DIR)/scale_land_grid_cartesC.o	\
	$(OBJ_DIR)/scale_land_grid_cartesC_index.o	\
	$(OBJ_DIR)/scale_land_grid_cartesC_real.o	\
	$(OBJ_DIR)/scale_land_grid_icoA.o	\
	$(OBJ_DIR)/scale_land_grid_icoA_index.o	\
	$(OBJ_DIR)/scale_land_dyn_bucket.o	\
	$(OBJ_DIR)/scale_land_phy_matsiro.o	\
	$(OBJ_DIR)/scale_land_phy_snow_diagnos.o	\
	$(OBJ_DIR)/scale_land_phy_snow_ky90.o	\
	$(OBJ_DIR)/scale_urban_grid_cartesC.o	\
	$(OBJ_DIR)/scale_urban_grid_cartesC_index.o	\
	$(OBJ_DIR)/scale_urban_grid_cartesC_real.o	\
	$(OBJ_DIR)/scale_urban_grid_icoA.o	\
	$(OBJ_DIR)/scale_urban_grid_icoA_index.o	\
	$(OBJ_DIR)/scale_urban_dyn_kusaka01.o	\
	$(OBJ_DIR)/scale_bulkflux.o	\
	$(OBJ_DIR)/scale_landuse.o	\
	$(OBJ_DIR)/scale_cpl_sfc_index.o	\
	$(OBJ_DIR)/scale_cpl_phy_sfc_fixed_temp.o	\
	$(OBJ_DIR)/scale_cpl_phy_sfc_skin.o

OBJS_FOR_DYN = \
	$(OBJ_DIR)/scale_atmos_dyn.o	\
	$(OBJ_DIR)/scale_atmos_dyn_common.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux_cd2.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux_cd4.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux_cd6.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux_ud1.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux_ud3.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux_ud3Koren1993.o	\
	$(OBJ_DIR)/scale_atmos_dyn_fvm_flux_ud5.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_large.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_large_euler.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_large_rk3.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_short.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_short_rk3.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_short_rk4.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_tracer.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_tracer_euler.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tinteg_tracer_rk3.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_large.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_large_fvm_heve.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_short.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_short_fvm_heve.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_short_fvm_hevi.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_short_fvm_hivi.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_tracer.o	\
	$(OBJ_DIR)/scale_atmos_dyn_tstep_tracer_fvm_heve.o

OBJS_FOR_SDM = \
	$(OBJ_DIR)/scale_precision.o	\
	$(OBJ_DIR)/scale_prof.o	\
	$(OBJ_DIR)/scale_comm_cartesC.o	\
	$(OBJ_DIR)/scale_atmos_thermodyn.o	\
	$(OBJ_DIR)/scale_atmos_saturation.o	\
	$(OBJ_DIR)/scale_atmos_phy_mp_common.o	\
	$(OBJ_DIR)/scale_tracer.o	\
	$(OBJ_DIR)/scale_tracer_sdm.o	\
	$(OBJ_DIR)/scale_io.o	\
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

$(LIBDIR)/$(LIBNAME): $(LIBNAME)
	install $< $@

$(LIBNAME): $(OBJS)
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
	rsync -th $(SDMDIR)/scale_atmos_phy_mp_sdm.F90 $(SCALELIBDIR)/src/atmosphere/physics/microphysics/
	@echo

restoresdm:
	@echo;echo "Restore SDM-related file"
	rsync -th $(SDMDIR)/scale_atmos_phy_mp_sdm_dummy.F90 $(SCALELIBDIR)/src/atmosphere/physics/microphysics/scale_atmos_phy_mp_sdm.F90
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
	rm -f $(MODDIR)/*.mod

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

file/scale_file.F90 : scale_file.F90.erb

file/scale_file_history.F90 : scale_file_history.F90.erb

scale_atmos_dyn_fvm_flux_%.F90 : scale_atmos_dyn_fvm_flux_udcd.F90.erb
	fname=$@ erb $< > $@

$(OBJ_DIR)/%.o : %.F90
	$(FC) $(FFLAGS)     $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(LIBVERSION)\" -I$(SCALELIBDIR)/include -I$(OBJ_DIR) $(SCALE_PAPI_INCLUDE) $(SCALE_NETCDF_INCLUDE) $(CONTRIB_INCLUDE) $(MODDIROPT) $(OBJ_DIR) -o $@ -c $<
$(OBJS_FOR_DYN) :
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(LIBVERSION)\" -I$(SCALELIBDIR)/include -I$(OBJ_DIR) $(SCALE_PAPI_INCLUDE) $(SCALE_NETCDF_INCLUDE) $(CONTRIB_INCLUDE) $(MODDIROPT) $(OBJ_DIR) -o $@ -c $<
$(OBJ_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -I$(SCALELIBDIR)/include $(SCALE_NETCDF_INCLUDE) -o $@ -c $<

.PHONY : clean distclean allclean depend modules dependency

include $(TOPDIR)/utils/make/Make_environments

ifneq ($(SCALE_IGNORE_SRCDEP),T)
-include depend
endif
# DO NOT DELETE
