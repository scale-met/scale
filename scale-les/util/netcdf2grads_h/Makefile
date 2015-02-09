##----------------------------------------------------------

TOPDIR      = $(abspath ../../..)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = ../../../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

LIBS=$(NETCDF_LIBS)
INCLUDES=$(NETCDF_INCLUDE)

OBJECT=net2g
MODS=prg_netcdf2grads_h.o

##----------------------------------------------------------
## make rules

all: main

main: prg_netcdf2grads_h.o
	$(FC) $(FFLAGS) -o $(OBJECT) $(MODS) $(LIBS)

#f90.o:
prg_netcdf2grads_h.o: prg_netcdf2grads_h.F90
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	rm -f *.o *.mod *.lst $(OBJECT)
##----------------------------------------------------------
