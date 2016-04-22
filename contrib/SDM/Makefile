#################################################################
#
# Makefile for Super-Droplet Method
#
#################################################################

TOPDIR      = ../..

SYSDEP_DIR = $(TOPDIR)/sysdep
include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/contrib/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

LIBNAME = libsdm.a 

SCALESRCDIR = $(SCALELIBDIR)/src

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

all: makedir $(LIBDIR)/$(LIBNAME) makemodules

makedir:
	@echo;echo "Entering contrib/SDM..."
	mkdir -p $(LIBDIR)

$(LIBDIR)/$(LIBNAME):	$(LIBNAME)
	mkdir -p $(LIBDIR)
	install $(LIBNAME) $(LIBDIR)

$(LIBNAME): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

makemodules: $(OBJS)

sdm_condensation_water.o: sdm_condensation_water.f90 \
	$(DEPENDLIB) \
	sdm_coordtrans.o \
	sdm_common.o \
	sdm_fluidconv.o \
	$(SCALESRCDIR)/scale_atmos_sub_saturation.o \
	$(SCALESRCDIR)/scale_tracer.o \
	$(SCALESRCDIR)/scale_const.o \
	$(SCALESRCDIR)/scale_grid_index.o \
	$(SCALESRCDIR)/scale_precision.o
	$(FC) $(FFLAGS) $(FFLAGS_SDM1) $(FFLAGS_SDM_MONIT) -I$(OBJ_DIR) -I$(SCALESRCDIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_sd2fluid.o: sdm_sd2fluid.f90 \
	$(DEPENDLIB) \
	sdm_common.o \
	sdm_coordtrans.o \
	$(SCALESRCDIR)/scale_const.o \
	$(SCALESRCDIR)/scale_precision.o
	$(FC) $(FFLAGS) $(FFLAGS_SDM2) $(FFLAGS_SDM_MONIT) -I$(OBJ_DIR) -I$(SCALESRCDIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_coalescence.o: sdm_coalescence.f90 \
	$(DEPENDLIB) \
	sdm_common.o \
	sdm_coordtrans.o \
	sdm_idutil.o \
	$(SCALESRCDIR)/scale_grid_index.o \
	$(SCALESRCDIR)/scale_const.o \
	$(SCALESRCDIR)/scale_precision.o \
	$(MTDIR)/rng_uniform_mt.mod \
	$(GAGETRYDIR)/gadg_algorithm.mod
	$(FC) $(FFLAGS) $(FFLAGS_SDM3) $(FFLAGS_SDM_MONIT) -I$(OBJ_DIR) -I$(SCALESRCDIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_idutil.o: sdm_idutil.f90 \
	$(DEPENDLIB) \
	sdm_common.o \
	$(SCALESRCDIR)/scale_grid_index.o \
	$(SCALESRCDIR)/scale_precision.o \
	$(GAGETRYDIR)/gadg_algorithm.mod
	$(FC) $(FFLAGS) $(FFLAGS_SDM_MONIT) -I$(OBJ_DIR) -I$(SCALESRCDIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_motion.o: sdm_motion.f90 \
	$(DEPENDLIB) \
	sdm_common.o \
	sdm_coordtrans.o \
	$(SCALESRCDIR)/scale_const.o \
	$(SCALESRCDIR)/scale_grid_index.o \
	$(SCALESRCDIR)/scale_grid_cartesian.o \
	$(SCALESRCDIR)/scale_precision.o
	$(FC) $(FFLAGS) $(FFLAGS_SDM4) $(FFLAGS_SDM_MONIT) -I$(OBJ_DIR) -I$(SCALESRCDIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_boundary.o: sdm_boundary.f90 \
	$(DEPENDLIB) \
	sdm_coordtrans.o \
	sdm_common.o \
	$(SCALESRCDIR)/scale_process.o \
	$(SCALESRCDIR)/scale_stdio.o \
	$(SCALESRCDIR)/scale_precision.o \
	$(SCALESRCDIR)/scale_grid_index.o \
	$(SCALESRCDIR)/scale_grid_cartesian.o
	$(FC) $(FFLAGS) $(FFLAGS_SDM_MONIT) -I$(OBJ_DIR) -I$(SCALESRCDIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_fluidconv.o: sdm_fluidconv.f90 \
	$(DEPENDLIB) \
	$(SCALESRCDIR)/scale_grid_cartesian.o \
	$(SCALESRCDIR)/scale_grid_index.o \
	$(SCALESRCDIR)/scale_const.o \
	$(SCALESRCDIR)/scale_atmos_sub_thermodyn.o \
	$(SCALESRCDIR)/scale_tracer.o \
	$(SCALESRCDIR)/scale_precision.o
	$(FC) $(FFLAGS) $(FFLAGS_SDM_MONIT) -I$(OBJ_DIR) -I$(SCALESRCDIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_io.o: sdm_io.f90 \
	$(DEPENDLIB) \
	$(SCALESRCDIR)/scale_process.o \
	$(SCALESRCDIR)/scale_precision.o \
	$(SCALESRCDIR)/scale_stdio.o
	$(FC) $(FFLAGS) $(FFLAGS_SDM_MONIT) -I$(OBJ_DIR) -I$(SCALESRCDIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_coordtrans.o: sdm_coordtrans.f90 \
	$(DEPENDLIB) \
	sdm_common.o \
	$(SCALESRCDIR)/scale_stdio.o \
	$(SCALESRCDIR)/scale_process.o \
	$(SCALESRCDIR)/scale_grid_index.o \
	$(SCALESRCDIR)/scale_grid_real.o \
	$(SCALESRCDIR)/scale_grid_cartesian.o \
	$(SCALESRCDIR)/scale_precision.o
	$(FC) $(FFLAGS) $(FFLAGS_SDM_MONIT) -I$(OBJ_DIR) -I$(SCALESRCDIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_memmgr.o: sdm_memmgr.f90 \
	$(DEPENDLIB) \
	sdm_common.o \
	$(SCALESRCDIR)/scale_tracer.o \
	$(SCALESRCDIR)/scale_grid_index.o \
	$(SCALESRCDIR)/scale_gridtrans.o
	$(FC) $(FFLAGS) $(FFLAGS_SDM_MONIT) -I$(OBJ_DIR) -I$(SCALESRCDIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_numset.o: sdm_numset.f90 \
	$(DEPENDLIB) \
	sdm_common.o \
	$(SCALESRCDIR)/scale_precision.o \
	$(SCALESRCDIR)/scale_grid_index.o \
	$(SCALESRCDIR)/scale_topography.o \
	$(SCALESRCDIR)/scale_grid_real.o
	$(FC) $(FFLAGS) $(FFLAGS_SDM_MONIT) -I$(OBJ_DIR) -I$(SCALESRCDIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

sdm_common.o: sdm_common.f90 \
	$(DEPENDLIB) \
	$(SCALESRCDIR)/scale_precision.o \
	$(SCALESRCDIR)/scale_stdio.o \
	$(SCALESRCDIR)/scale_const.o \
	$(MTDIR)/rng_uniform_mt.mod
	$(FC) $(FFLAGS) $(FFLAGS_SDM_MONIT) -I$(OBJ_DIR) -I$(SCALESRCDIR) -I$(MODDIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MTDIR) -I$(GAGETRYDIR) -I$(SDMDIR) -o $@ -c $<

clean:
	rm -f $(MODS)
	rm -f $(OBJS)

distclean: clean
	rm -f $(LIBNAME)

allclean: distclean
	rm -f $(LIBDIR)/$(LIBNAME)
	rm -f *.lst
	cd $(MODDIR); rm -f $(MODS)

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.F90.o:
	$(FC) $(FFLAGS) -o $@ -c $<
.F95.o:
	$(FC) $(FFLAGS) -o $@ -c $<

%.o: %.mod

.PHONY : clean distclean allclean
