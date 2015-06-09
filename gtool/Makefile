################################################################################
#
# Makefile for main program
#
################################################################################

TOPDIR      = ..
SYSDEP_DIR  = ../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

LIBNAME = libgtool.a

DEPENDLIB = $(LIBDIR)/libdcutils.a

MODS = \
	gtool_file_h.o  \
	gtool_file.o    \
	gtool_history.o \
	gtool_file_f.o  \
	gtool_netcdf.o

all: makedir $(LIBDIR)/$(LIBNAME)

makedir:
	@echo;echo "Entering gtool..."
	mkdir -p $(LIBDIR)

$(LIBDIR)/libdcutils.a:
	$(MAKE) -C $(DCUTILSDIR)

$(LIBDIR)/$(LIBNAME): $(LIBNAME)
	install $(LIBNAME) $(LIBDIR)

$(LIBNAME): $(MODS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

gtool_file_h.o  : gtool_file_h.f90  $(DEPENDLIB)
gtool_file.o    : gtool_file.f90    $(DEPENDLIB) gtool_file_h.o
gtool_history.o : gtool_history.f90 $(DEPENDLIB) gtool_file_h.o gtool_file.o
gtool_file_f.o  : gtool_file_f.c
gtool_netcdf.o  : gtool_netcdf.c



allclean: distclean
	$(MAKE) -C $(DCUTILSDIR) distclean
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
	$(FC) $(FFLAGS) -I$(DCUTILSDIR) -o $@ -c $<
.c.o:
	$(CC) $(CFLAGS) $(NETCDF_INCLUDE) -o $@ -c $<

%.o: %.mod

.PHONY : clean distclean allclean
