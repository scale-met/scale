################################################################################
#
# Makefile for main program
#
################################################################################


TOPDIR = ..

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

MODS = \
	dc_calendar.o \
	dc_log.o


all: makelibdcutils

makelibdcutils: libdcutils.a
	mkdir -p $(LIBDIR)
	install libdcutils.a $(LIBDIR)

libdcutils.a: $(MODS)
	$(AR) $@ $?
	$(RANLIB) $@

dc_calendar.f90: dc_calendar.f90.erb
	erb $< > $@

dc_calendar.o: dc_calendar.f90 dc_log.o

clean:
	rm -f $(MODS) *.mod
distclean: clean
	rm -f libdcutils.a

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.f90.o:
	$(FC) $(FFLAGS) -o $@ -c $<
.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

%.o: %.mod

.PHONY : clean distclean
