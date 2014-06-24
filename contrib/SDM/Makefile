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

MODS	= \
	  m_sdm_memmgr.mod \
	  m_sdm_numset.mod \
	  m_sdm_common.mod

OBJS	= \
	  sdm_memmgr.o \
	  sdm_numset.o \
	  sdm_common.o

all:	$(LIBDIR)/$(LIBNAME) instmodules

$(LIBDIR)/$(LIBNAME):	$(LIBNAME)
	mkdir -p $(LIBDIR)
	install $(LIBNAME) $(LIBDIR)

$(LIBNAME): $(OBJS)
	$(AR) r $@ $?
	$(RANLIB) $@

instmodules:	makemodules
	mkdir -p $(MODDIR)
	install $(MODS) $(MODDIR)

makemodules: $(OBJS)

sdm_memmgr.o: sdm_memmgr.f90 \
	sdm_common.o \
	$(BUILD_DIR)/scale_grid_index.o \
	$(BUILD_DIR)/scale_gridtrans.o
	$(FC) $(FFLAGS) -I$(BUILD_DIR) -I$(MODDIR) -o $@ -c $<

sdm_numset.o: sdm_numset.f90 \
	sdm_common.o \
	$(BUILD_DIR)/scale_precision.o \
	$(BUILD_DIR)/scale_grid_index.o \
	$(BUILD_DIR)/scale_topography.o \
	$(BUILD_DIR)/scale_grid_real.o 
	$(FC) $(FFLAGS) -I$(BUILD_DIR) -I$(MODDIR) -o $@ -c $<

sdm_common.o: sdm_common.f90 \
	$(BUILD_DIR)/scale_precision.o \
	$(BUILD_DIR)/scale_stdio.o \
	$(BUILD_DIR)/scale_const.o \
	$(MODDIR)/rng_uniform_mt.mod
	$(FC) $(FFLAGS) -I$(BUILD_DIR) -I$(MODDIR) -o $@ -c $<

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
	$(FC) $(FFLAGS) -o $@ -c $<
.F95.o:
	$(FC) $(FFLAGS) -o $@ -c $<

%.o: %.mod

.PHONY : clean distclean allclean
