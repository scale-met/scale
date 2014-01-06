################################################################################
#
# Makefile for main program
#
################################################################################

TOPDIR      = ..
SYSDEP_DIR  = ../sysdep

include $(TOPDIR)/Mkinclude
include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)

LIBNAME = libgtool.a

MODS = \
	gtool_file.o \
	gtool_file_f.o \
	gtool_file_h.o \
	gtool_history.o \
	gtool_netcdf.o


all: $(LIBDIR)/$(LIBNAME)

$(LIBDIR)/$(LIBNAME): $(LIBNAME)
	mkdir -p $(LIBDIR)
	install $(LIBNAME) $(LIBDIR)

$(LIBNAME): $(DCUTILSDIR)/libdcutils.a $(MODS)
	$(AR) r $@ $(MODS)
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

clean:
	rm -f $(MODS) *.mod
distclean: clean
	rm -f $(LIBNAME)
allclean: distclean
	cd $(DCUTILSDIR); $(MAKE) distclean
	rm -f $(LIBDIR)/$(LIBNAME)

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.f90.o:
	$(FC) $(FFLAGS) -I$(DCUTILSDIR) -o $@ -c $<

%.o: %.mod

.PHONY : clean distclean allclean
