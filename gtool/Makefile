################################################################################
#
# Makefile for main program
#
################################################################################


include ../sysdep/Makedef.$(SCALE_SYS)

LIBDIR = ../lib

MODS = \
	gtool_file.o \
	gtool_file_f.o \
	gtool_file_h.o \
	gtool_history.o \
	gtool_netcdf.o


all: $(LIBDIR)/libgtool.a

$(LIBDIR)/libgtool.a: libgtool.a
	mkdir -p $(LIBDIR)
	install libgtool.a $(LIBDIR)

libgtool.a: $(LIBDIR)/libdcutils.a $(MODS)
	ar r $@ $?

gtool_file.f90: gtool_file.f90.erb
	erb $< > $@

gtool_history.f90: gtool_history.f90.erb
	erb $< > $@

gtool_netcdf.o: gtool_netcdf.c
	$(CC) $(CFLAGS) $(NETCDF_INCLUDE) -o $@ -c $<

gtool_file.o: gtool_file.f90 gtool_file_h.o dc_log.o
gtool_history.o: gtool_history.f90 dc_log.o gtool_file_h.o

clean:
	rm -f $(MODS) *.mod
distclean: clean
	rm -f libgtool.a

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.f90.o:
	$(FC) $(FFLAGS) -o $@ -c $<
.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

%.o: %.mod

.PHONY : clean distclean
