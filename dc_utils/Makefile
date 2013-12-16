################################################################################
#
# Makefile for main program
#
################################################################################

TOPDIR      = ..

include $(TOPDIR)/Mkinclude
include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)

LIBNAME = libdcutils.a

MODS = \
	dc_calendar.o \
	dc_log.o \
	dc_types.o \
	dc_test.o

all: $(LIBDIR)/$(LIBNAME)

$(LIBDIR)/$(LIBNAME): $(LIBNAME)
	mkdir -p $(LIBDIR)
	install $(LIBNAME) $(LIBDIR)

$(LIBNAME): $(MODS)
	$(AR) r $@ $?
	$(RANLIB) $@

dc_calendar.f90: dc_calendar.f90.erb
	erb $< > $@

dc_test.f90: dc_test.rb2f90
	ruby $< > $@

dc_calendar.o: dc_calendar.f90 dc_log.o dc_types.o

dc_types.o: dc_types.f90

dc_test.o: dc_test.f90 dc_types.o
	$(FC) $(FFLAGS_DEBUG) -o $@ -c $<

clean:
	rm -f $(MODS) *.mod
distclean: clean
	rm -f $(LIBNAME)
allclean: distclean
	rm -f $(LIBDIR)/$(LIBNAME)

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.f90.o:
	$(FC) $(FFLAGS) -o $@ -c $<

%.o: %.mod

.PHONY : clean distclean allclean
