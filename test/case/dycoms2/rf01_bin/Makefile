################################################################################
#
# Makefile for each test program
#
################################################################################

TOPDIR      = ../../../..
BUILD_DIR  := $(shell pwd)/.libs
CODE_DIR := $(shell pwd)/code
INDEX_DIR  := $(shell pwd)
SYSDEP_DIR := $(shell pwd)/$(TOPDIR)/sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

# for build
index      = 200x8x8
precision  = real8
dynamics   = heve
tracer     = binw
sfcflx     = const
turbulence = smg
radiation  = dycoms2
postfix    = 

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 6



all: buildlocal build 

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common



buildlocal: makedir $(BUILD_DIR)/mod_atmos_phy_rd_dycoms2.f90 $(BUILD_DIR)/mod_user.f90

$(BUILD_DIR)/mod_atmos_phy_rd_dycoms2.f90: $(CODE_DIR)/mod_atmos_phy_rd_dycoms2.f90
$(BUILD_DIR)/mod_user.f90: $(CODE_DIR)/mod_user.f90
	cp -f $< $@

analy: 
	$(FC) $(FFLAGS) -L$(NETCDF_LIBS) -lnetcdff -lnetcdf -I$(NETCDF_INCLUDE) -o $(BUILD_DIR)/profile_maker $(CODE_DIR)/make_profile_cmp_nc_each.f90
