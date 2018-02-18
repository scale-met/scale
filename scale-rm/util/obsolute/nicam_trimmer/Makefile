##----------------------------------------------------------

TOPDIR      = $(abspath ../../..)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = $(TOPDIR)/sysdep
DCUTILS     = $(TOPDIR)/dc_utils

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

#FC=ifort
LIBS=$(SCALE_NETCDF_LIBS) -L$(TOPDIR)/lib -lscale  -L$(DCUTILS) -ldcutils
INCLUDES=-I$(TOPDIR)/include -I$(DCUTILS) $(SCALE_NETCDF_INCLUDE)

OBJECT1=nicam_trimmer
OBJECT2=nicam_trimmer_jw06
MODS1=mod_netcdf.o prg_nicam_trimmer.o
MODS2=mod_netcdf.o prg_nicam_trimmer_jw06.o

##----------------------------------------------------------
## make rules

all: main jw06

main: prg_nicam_trimmer.o
	$(FC) $(FFLAGS) -o $(OBJECT1) $(MODS1) $(LIBS)

jw06: prg_nicam_trimmer_jw06.o
	$(FC) $(FFLAGS) -o $(OBJECT2) $(MODS2) $(LIBS)

#f90.o:
mod_netcdf.o: src/mod_netcdf.f90
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
prg_nicam_trimmer.o: src/prg_nicam_trimmer.f90 mod_netcdf.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
prg_nicam_trimmer_jw06.o: src/prg_nicam_trimmer_jw06.f90 mod_netcdf.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	rm -f *.o *.mod $(OBJECT)
##----------------------------------------------------------
