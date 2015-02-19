##----------------------------------------------------------

TOPDIR      = $(abspath ../../..)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = ../../../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)

LIBS=$(NETCDF_LIBS)
INCLUDES=$(NETCDF_INCLUDE)
#FFLAGS=$(FFLAGS_FAST)

OBJECT=popsca
MODS=prg_popsca_h.o

##----------------------------------------------------------
## make rules

all: main

main: prg_popsca_h.o
	$(FC) $(FFLAGS) -o $(OBJECT) $(MODS) $(LIBS)

#f90.o:
prg_popsca_h.o: prg_popsca_h.F90
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	rm -f *.o *.mod *.lst $(OBJECT)
##----------------------------------------------------------
