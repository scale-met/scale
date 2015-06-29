##----------------------------------------------------------

TOPDIR      = $(abspath ../../..)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = ../../../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

LIBS=$(NETCDF_LIBS)
INCLUDES=$(NETCDF_INCLUDE) -Iinc

OBJECT=net2g

ifeq ($(NOMPI),T)
  FC=$(SFC)
  FFLAGS += -UMPIUSE
else
  FFLAGS += -DMPIUSE
endif

MODS=mod_net2g_vars.o     \
     mod_net2g_error.o    \
     mod_net2g_setup.o    \
     mod_net2g_anal.o     \
     mod_net2g_calender.o \
     mod_net2g_comm.o     \
     mod_net2g_netcdf.o   \
     prg_netcdf2grads_h.o

##----------------------------------------------------------
## make rules

all: main

main: prg_netcdf2grads_h.o
	$(FC) $(FFLAGS) -o $(OBJECT) $(MODS) $(LIBS)
	@mv *.mod inc/; mv *.o inc/;

#f90.o:
mod_net2g_vars.o: src/mod_net2g_vars.F90
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_net2g_error.o: src/mod_net2g_error.F90 mod_net2g_vars.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_net2g_setup.o: src/mod_net2g_setup.F90 mod_net2g_vars.o mod_net2g_error.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_net2g_anal.o: src/mod_net2g_anal.F90 mod_net2g_vars.o mod_net2g_error.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_net2g_calender.o: src/mod_net2g_calender.F90 mod_net2g_vars.o mod_net2g_error.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_net2g_comm.o: src/mod_net2g_comm.F90 mod_net2g_vars.o mod_net2g_error.o mod_net2g_setup.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_net2g_netcdf.o: src/mod_net2g_netcdf.F90 mod_net2g_vars.o mod_net2g_error.o mod_net2g_anal.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
prg_netcdf2grads_h.o: src/prg_netcdf2grads_h.F90 mod_net2g_vars.o mod_net2g_error.o \
        mod_net2g_setup.o mod_net2g_anal.o mod_net2g_calender.o mod_net2g_comm.o mod_net2g_netcdf.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	rm -f inc/*.o inc/*.mod *.lst $(OBJECT)
##----------------------------------------------------------
