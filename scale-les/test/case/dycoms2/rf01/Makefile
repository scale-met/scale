################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

# for build
index      = 200x8x8
precision  = real8
dynamics   = heve
tracer     = sn13
sfcflx     = const
turbulence = smg
radiation  = dycoms2

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 6

# srcfiles
CODE_DIR = $(PWD)/code
ORG_SRCS = mod_user.f90


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

analy:
	$(FC) $(FFLAGS_FAST) -o .libs/profile_maker \
	$(CODE_DIR)/make_profile_cmp_nc_each.f90 $(NETCDF_INCLUDE) $(NETCDF_LIBS)
