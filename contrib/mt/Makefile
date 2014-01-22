#################################################################
#
# Makefile for Mersenne Twister
#
#################################################################

TOPDIR      = ../..

include $(TOPDIR)/Mkinclude
include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)

LIBNAME = libmt.a 

MODS	= \
	  rng_uniform_mt.o

all:	$(LIBDIR)/$(LIBNAME)

$(LIBDIR)/$(LIBNAME):	$(LIBNAME)
	mkdir -p $(LIBDIR)
	install $(LIBNAME) $(LIBDIR)

$(LIBNAME): $(MODS)
	$(AR) r $@ $?
	$(RANLIB) $@

rng_uniform_mt.o: rng_uniform_mt.F90
	$(FC) $(FFLAGS) -o $@ -c $<

clean:
	rm -f $(MODS) *.mod

distclean: clean
	rm -f $(LIBNAME)

allclean: distclean
	rm -f $(LIBDIR)/$(LIBNAME)

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.F90.o:
	$(FC) $(FFLAGS) -o $@ -c $<
.F95.o:
	$(FC) $(FFLAGS) -o $@ -c $<

%.o: %.mod

.PHONY : clean distclean allclean
