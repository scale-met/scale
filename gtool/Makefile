################################################################################
#
# Makefile for main program
#
################################################################################


TOPDIR = ..

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

MODS = \
	gtool_file.o \
	gtool_file_f.o \
	gtool_file_h.o \
	gtool_history.o \
	gtool_netcdf.o


all: makelibgtool

makelibgtool: libgtool.a
	mkdir -p $(LIBDIR)
	install libgtool.a $(LIBDIR)

libgtool.a: $(DCUTILSDIR)/libdcutils.a $(MODS)
	$(AR) $@ $(MODS)
	$(RANLIB) $@

$(DCUTILSDIR)/libdcutils.a:
	cd $(DCUTILSDIR); $(MAKE)

gtool_file.f90: gtool_file.f90.erb
	erb $< > $@

gtool_history.f90: gtool_history.f90.erb
	erb $< > $@

gtool_netcdf.o: gtool_netcdf.c
	$(CC) $(CFLAGS) $(NETCDF_INCLUDE) -o $@ -c $<

gtool_file.o: gtool_file.f90 gtool_file_h.o $(LIBDIR)/libdcutils.a
gtool_history.o: gtool_history.f90 gtool_file_h.o $(LIBDIR)/libdcutils.a

allclean: distclean
	cd $(DCUTILSDIR); $(MAKE) distclean
clean:
	rm -f $(MODS) *.mod
distclean: clean
	rm -f libgtool.a

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.f90.o:
	$(FC) $(FFLAGS) -I$(DCUTILSDIR) -o $@ -c $<

%.o: %.mod

.PHONY : clean distclean allclean
