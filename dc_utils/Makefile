################################################################################
#
# Makefile for main program
#
################################################################################

TOPDIR      = ..
SYSDEP_DIR  = ../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

LIBNAME = libdcutils.a

MODS = \
	dc_log.o      \
	dc_types.o    \
	dc_calendar.o \
	dc_test.o



all: makedir $(LIBDIR)/$(LIBNAME)

makedir:
	@echo;echo "Entering dcutils..."
	mkdir -p $(LIBDIR)

$(LIBDIR)/$(LIBNAME): $(LIBNAME)
	mkdir -p $(LIBDIR)
	install $(LIBNAME) $(LIBDIR)

$(LIBNAME): $(MODS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

dc_types.o    : dc_types.f90
dc_log.o      : dc_log.f90
dc_calendar.o : dc_calendar.f90 dc_types.o dc_log.o
dc_test.o     : dc_test.f90     dc_types.o
	$(FC) $(FFLAGS_DEBUG) -o $@ -c $<

dc_test.f90: dc_test.rb2f90
	ruby $< > $@



allclean: distclean
	rm -f $(LIBDIR)/$(LIBNAME)

distclean: clean
	rm -f $(LIBNAME)

clean:
	rm -f $(MODS) *.mod



.SUFFIXES:
.SUFFIXES: .o .f90 .c .erb .mod

%.f90 : %.f90.erb
	erb $< > $@

.f90.o:
	$(FC) $(FFLAGS) -o $@ -c $<

%.o: %.mod

.PHONY : clean distclean allclean
