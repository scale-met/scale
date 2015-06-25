#################################################################
#
# Makefile for Super-Droplet Method
#
#################################################################

TOPDIR      = ../..

SYSDEP_DIR = $(TOPDIR)/sysdep
include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

LIBNAME = libsdm.a 

DEPENDLIB = $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a

MODS	= \
	  m_sdm_condensation_water.mod \
	  m_sdm_sd2fluid.mod \
	  m_sdm_coalescence.mod \
	  m_sdm_idutil.mod \
	  m_sdm_motion.mod \
	  m_sdm_boundary.mod \
	  m_sdm_fluidconv.mod \
	  m_sdm_io.mod \
	  m_sdm_coordtrans.mod \
	  m_sdm_memmgr.mod \
	  m_sdm_numset.mod \
	  m_sdm_common.mod

OBJS	= \
	  sdm_condensation_water.o \
	  sdm_sd2fluid.o \
	  sdm_coalescence.o \
	  sdm_idutil.o \
	  sdm_motion.o \
	  sdm_boundary.o \
	  sdm_fluidconv.o \
	  sdm_io.o \
	  sdm_coordtrans.o \
	  sdm_memmgr.o \
	  sdm_numset.o \
	  sdm_common.o

all:	$(LIBDIR)/$(LIBNAME) makemodules

$(LIBDIR)/$(LIBNAME):	$(LIBNAME)
	mkdir -p $(LIBDIR)
	install $(LIBNAME) $(LIBDIR)

$(LIBNAME): $(OBJS)
	$(AR) r $@ $?
	$(RANLIB) $@

makemodules: $(OBJS)

sdm_condensation_water.o: sdm_condensation_water.f90 \
	$(DEPENDLIB) \
	sdm_coordtrans.o \
	sdm_common.o \
	sdm_fluidconv.o \
        $(BUILD_DIR)/scale_atmos_sub_saturation.o \
	$(BUILD_DIR)/scale_tracer.o \
	$(BUILD_DIR)/scale_const.o \
	$(BUILD_DIR)/scale_grid_index.o \
	$(BUILD_DIR)/scale_precision.o
	$(FC) $(FFLAGS_SDM) -Kprefetch_indirect,simd=2 $(FFLAGS_SDM_MONIT) -I$(BUILD_DIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_sd2fluid.o: sdm_sd2fluid.f90 \
	$(DEPENDLIB) \
	sdm_common.o \
	sdm_coordtrans.o \
	$(BUILD_DIR)/scale_const.o \
	$(BUILD_DIR)/scale_precision.o
	$(FC) $(FFLAGS_SDM) $(FFLAGS_SDM_MONIT) -I$(BUILD_DIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_coalescence.o: sdm_coalescence.f90 \
	$(DEPENDLIB) \
	sdm_common.o \
	sdm_coordtrans.o \
	sdm_idutil.o \
	$(BUILD_DIR)/scale_grid_index.o \
	$(BUILD_DIR)/scale_const.o \
	$(BUILD_DIR)/scale_precision.o \
	$(MTDIR)/rng_uniform_mt.mod \
	$(GAGETRYDIR)/gadg_algorithm.mod
	$(FC) $(FFLAGS_SDM) -Ksimd=2 $(FFLAGS_SDM_MONIT) -I$(BUILD_DIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_idutil.o: sdm_idutil.f90 \
	$(DEPENDLIB) \
	sdm_common.o \
	$(BUILD_DIR)/scale_grid_index.o \
	$(BUILD_DIR)/scale_precision.o \
	$(GAGETRYDIR)/gadg_algorithm.mod
	$(FC) $(FFLAGS_SDM) $(FFLAGS_SDM_MONIT) -I$(BUILD_DIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_motion.o: sdm_motion.f90 \
	$(DEPENDLIB) \
	sdm_common.o \
	sdm_coordtrans.o \
	$(BUILD_DIR)/scale_const.o \
	$(BUILD_DIR)/scale_grid_index.o \
	$(BUILD_DIR)/scale_grid_cartesian.o \
	$(BUILD_DIR)/scale_precision.o
	$(FC) $(FFLAGS_SDM) -Ksimd=2 $(FFLAGS_SDM_MONIT) -I$(BUILD_DIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_boundary.o: sdm_boundary.f90 \
	$(DEPENDLIB) \
	sdm_coordtrans.o \
	sdm_common.o \
	$(BUILD_DIR)/scale_process.o \
	$(BUILD_DIR)/scale_stdio.o \
	$(BUILD_DIR)/scale_precision.o \
	$(BUILD_DIR)/scale_grid_index.o \
	$(BUILD_DIR)/scale_grid_cartesian.o
	$(FC) $(FFLAGS_SDM) $(FFLAGS_SDM_MONIT) -I$(BUILD_DIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_fluidconv.o: sdm_fluidconv.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/scale_grid_cartesian.o \
	$(BUILD_DIR)/scale_grid_index.o \
	$(BUILD_DIR)/scale_const.o \
	$(BUILD_DIR)/scale_atmos_sub_thermodyn.o \
	$(BUILD_DIR)/scale_tracer.o \
	$(BUILD_DIR)/scale_precision.o
	$(FC) $(FFLAGS_SDM) $(FFLAGS_SDM_MONIT) -I$(BUILD_DIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_io.o: sdm_io.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/scale_process.o \
	$(BUILD_DIR)/scale_precision.o \
	$(BUILD_DIR)/scale_stdio.o
	$(FC) $(FFLAGS_SDM) $(FFLAGS_SDM_MONIT) -I$(BUILD_DIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_coordtrans.o: sdm_coordtrans.f90 \
	$(DEPENDLIB) \
	sdm_common.o \
	$(BUILD_DIR)/scale_stdio.o \
	$(BUILD_DIR)/scale_process.o \
	$(BUILD_DIR)/scale_grid_index.o \
	$(BUILD_DIR)/scale_grid_real.o \
	$(BUILD_DIR)/scale_grid_cartesian.o \
	$(BUILD_DIR)/scale_precision.o
	$(FC) $(FFLAGS_SDM) $(FFLAGS_SDM_MONIT) -I$(BUILD_DIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_memmgr.o: sdm_memmgr.f90 \
	$(DEPENDLIB) \
	sdm_common.o \
	$(BUILD_DIR)/scale_tracer.o \
	$(BUILD_DIR)/scale_grid_index.o \
	$(BUILD_DIR)/scale_gridtrans.o
	$(FC) $(FFLAGS_SDM) $(FFLAGS_SDM_MONIT) -I$(BUILD_DIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_numset.o: sdm_numset.f90 \
	$(DEPENDLIB) \
	sdm_common.o \
	$(BUILD_DIR)/scale_precision.o \
	$(BUILD_DIR)/scale_grid_index.o \
	$(BUILD_DIR)/scale_topography.o \
	$(BUILD_DIR)/scale_grid_real.o 
	$(FC) $(FFLAGS_SDM) $(FFLAGS_SDM_MONIT) -I$(BUILD_DIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_common.o: sdm_common.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/scale_precision.o \
	$(BUILD_DIR)/scale_stdio.o \
	$(BUILD_DIR)/scale_const.o \
	$(MTDIR)/rng_uniform_mt.mod
	$(FC) $(FFLAGS_SDM) $(FFLAGS_SDM_MONIT) -I$(BUILD_DIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

clean:
	rm -f $(MODS)
	rm -f $(OBJS)

distclean: clean
	rm -f $(LIBNAME)

allclean: distclean
	rm -f $(LIBDIR)/$(LIBNAME)
	cd $(MODDIR); rm -f $(MODS)

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.F90.o:
	$(FC) $(FFLAGS_SDM) -o $@ -c $<
.F95.o:
	$(FC) $(FFLAGS_SDM) -o $@ -c $<

%.o: %.mod

.PHONY : clean distclean allclean
