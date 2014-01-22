#################################################################
# Makefile for Gadgetry
#################################################################
#SHELL=/bin/sh

TOPDIR      = ../..

include $(TOPDIR)/Mkinclude
include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)

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

all:	$(LIBDIR)/$(LIBNAME)

$(LIBDIR)/$(LIBNAME):	$(LIBNAME)
	mkdir -p $(LIBDIR)
	install $(LIBNAME) $(LIBDIR)

$(LIBNAME):	$(OBJECTS)
	$(AR) r $@ $?
	$(RANLIB) $@

include depend.mk

clean:
	rm -f $(OBJECTS)
	rm -f $(MODULES)
	rm -f $(LIBNAME)

distclean: clean
	rm -f $(LIBNAME)

allclean:
	rm -f $(OBJECTS)
	rm -f $(MODULES)
	rm -f $(LIBNAME)
	rm -f $(LIBDIR)/$(LIBNAME)

include depend.mk
