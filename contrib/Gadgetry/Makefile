#################################################################
# Makefile for Gadgetry
#################################################################
#SHELL=/bin/sh

TOPDIR      = ../..

SYSDEP_DIR = $(TOPDIR)/sysdep
include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

LIBNAME = libgagetry.a 

.PHONY:	 r dist depend
include files.mk

FPX3    = perl tools/fpx3 -I.

.PHONY:	all
.SUFFIXES: .o .F90 .F90_fpx3 .inc_fpx3 

.F90.o:
	$(FC) $(FFLAGS) -c $< -o $@

.F90_fpx3.F90:
	$(FPX3) $< > $@

.inc_fpx3.inc:
	$(FPX3) $< > $@

all:	$(LIBDIR)/$(LIBNAME) makemodules

$(LIBDIR)/$(LIBNAME):	$(LIBNAME)
	mkdir -p $(LIBDIR)
	install $(LIBNAME) $(LIBDIR)

$(LIBNAME):	$(OBJECTS)
	$(AR) r $@ $?
	$(RANLIB) $@

makemodules: $(OBJECTS)

include depend.mk

clean:
	rm -f $(OBJECTS)
	rm -f $(MODULES)

distclean: clean
	rm -f $(LIBNAME)

allclean:
	rm -f $(OBJECTS)
	rm -f $(MODULES)
	rm -f $(LIBNAME)
	rm -f $(LIBDIR)/$(LIBNAME)
	cd $(MODDIR); rm -f $(MODULES)
