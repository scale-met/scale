#################################################################
#
# Makefile for Mersenne Twister
#
#################################################################

TOPDIR      = ../..

SYSDEP_DIR = $(TOPDIR)/sysdep
include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/contrib/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

LIBNAME = libmt.a 

MODS	= \
	  rng_uniform_mt.mod

OBJS	= \
	  rng_uniform_mt.o

all: makedir $(LIBDIR)/$(LIBNAME) makemodules

makedir:
	@echo;echo "Entering contrib/mt..."
	mkdir -p $(LIBDIR)

$(LIBDIR)/$(LIBNAME): $(LIBNAME)
	mkdir -p $(LIBDIR)
	install $(LIBNAME) $(LIBDIR)

$(LIBNAME): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

makemodules: $(OBJS)

rng_uniform_mt.o: rng_uniform_mt.F90
	$(FC) $(FFLAGS) -o $@ -c $<

clean:
	rm -f $(MODS)
	rm -f $(OBJS)

distclean: clean
	rm -f $(LIBNAME)

allclean: distclean
	rm -f $(LIBDIR)/$(LIBNAME)
	cd $(MODDIR); rm -f $(MODS); rm -f *.lst

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.F90.o:
	$(FC) $(FFLAGS) -o $@ -c $<
.F95.o:
	$(FC) $(FFLAGS) -o $@ -c $<

%.o: %.mod

.PHONY : clean distclean allclean
