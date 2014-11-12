##----------------------------------------------------------

TOPDIR      = $(abspath ../../..)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = ../../../sysdep
SCALELIB    = ../../../scalelib/src
GTOOLLIB    = ../../../gtool
DCUTILS     = ../../../dc_utils

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)

#FC=ifort
LIBS=$(NETCDF_LIBS) -L$(SCALELIB) -lscale  -L$(DCUTILS) -ldcutils -L$(GTOOLLIB) -lgtool
INCLUDES=-I$(SCALELIB) -I$(DCUTILS) -I$(GTOOLLIB) $(NETCDF_INCLUDE)

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
mod_netcdf.o: mod_netcdf.f90
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
prg_nicam_trimmer.o: prg_nicam_trimmer.f90 mod_netcdf.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
prg_nicam_trimmer_jw06.o: prg_nicam_trimmer_jw06.f90 mod_netcdf.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	rm -f *.o *.mod $(OBJECT)
##----------------------------------------------------------
