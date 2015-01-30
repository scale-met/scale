##----------------------------------------------------------

TOPDIR      = $(abspath ../../..)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = ../../../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)

#FC=ifort
LIBS=$(NETCDF_LIBS)
INCLUDES=$(NETCDF_INCLUDE)

OBJECT=net2g
MODS=prg_net2g.o

##----------------------------------------------------------
## make rules

all: main

main: prg_net2g.o
	$(FC) $(FFLAGS) -o $(OBJECT) $(MODS) $(LIBS)

#f90.o:
prg_net2g.o: prg_net2g.F90
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	rm -f *.o *.mod $(OBJECT)
##----------------------------------------------------------
