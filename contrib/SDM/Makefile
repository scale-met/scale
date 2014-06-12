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
	  sdm_common.mod
OBJS	= \
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

sdm_common.o: sdm_common.f90 \
	$(BUILD_DIR)/scale_precision.mod \
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
