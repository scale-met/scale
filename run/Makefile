################################################################################
#
# Makefile for main program
#
################################################################################

SRCDIR = ../src
include $(SRCDIR)/sysdep/Makedef.$(SCALE_SYS)

LIBDIR = ../lib

all: check

makegtool:
	cd $(SRCDIR)/gtool; make
	sleep 1

check: check.o
	$(LD) $(LFLAGS) -o $@ $< -L$(LIBDIR) -lgtool $(NETCDF_LIBS)

check.o: check.f90 $(LIBDIR)/libgtool.a

$(LIBDIR)/libgtool.a: makegtool


clean:
	rm -f check.o

distclean: clean
	rm -f check

.SUFFIXES:
.SUFFIXES: .o .f90

.f90.o:
	$(FC) $(FFLAGS) -I$(SRCDIR)/include -I$(SRCDIR)/gtool -o $@ -c $<

.PHONY : clean distclean
