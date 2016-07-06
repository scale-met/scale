##----------------------------------------------------------

TOPDIR      = $(abspath ../../..)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = ../../../sysdep
SCALELIB    = ../../../lib
SCALEINC    = ../../../scalelib/src
GTOOLLIB    = ../../../gtool
DCUTILS     = ../../../dc_utils

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

LIBS=$(NETCDF_LIBS) -L$(SCALELIB) -lscale -lgtool -ldcutils
INCLUDES=$(NETCDF_INCLUDE) -Iinc -I$(SCALEINC) -I$(DCUTILS) -I$(GTOOLLIB)

OBJECT=net2g
MODS=mod_net2g_vars.o     \
     mod_net2g_error.o    \
     mod_net2g_setup.o    \
     mod_net2g_io.o       \
     mod_net2g_anal.o     \
     mod_net2g_calender.o \
     mod_net2g_comm.o     \
     mod_net2g_netcdf.o   \
     mod_netcdf2grads_h.o \
     prg_bulk_net2g.o

##----------------------------------------------------------
## make rules

all: main

main: prg_bulk_net2g.o
	$(FC) $(FFLAGS) -o $(OBJECT) $(MODS) $(LIBS)
	@mv *.mod inc/; mv *.o inc/;

#f90.o:
mod_net2g_vars.o: src/mod_net2g_vars.F90
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_net2g_error.o: src/mod_net2g_error2.F90 mod_net2g_vars.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_net2g_setup.o: src/mod_net2g_setup.F90 mod_net2g_vars.o mod_net2g_error.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_net2g_io.o: src/mod_net2g_io.F90 mod_net2g_vars.o mod_net2g_error.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_net2g_anal.o: src/mod_net2g_anal.F90 mod_net2g_vars.o mod_net2g_error.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_net2g_calender.o: src/mod_net2g_calender.F90 mod_net2g_vars.o mod_net2g_error.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_net2g_comm.o: src/mod_net2g_comm2.F90 mod_net2g_vars.o mod_net2g_error.o mod_net2g_setup.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_net2g_netcdf.o: src/mod_net2g_netcdf.F90 mod_net2g_vars.o mod_net2g_error.o mod_net2g_anal.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_netcdf2grads_h.o: src/mod_netcdf2grads_h.F90 mod_net2g_vars.o mod_net2g_error.o mod_net2g_setup.o \
        mod_net2g_io.o mod_net2g_anal.o mod_net2g_calender.o mod_net2g_comm.o mod_net2g_netcdf.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
prg_bulk_net2g.o: src/prg_bulk_net2g.f90 mod_net2g_vars.o mod_net2g_error.o mod_net2g_setup.o \
        mod_net2g_io.o mod_net2g_anal.o mod_net2g_calender.o mod_net2g_comm.o mod_net2g_netcdf.o mod_netcdf2grads_h.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	rm -f *.o *.mod inc/*.o inc/*.mod *.lst $(OBJECT)
##----------------------------------------------------------
