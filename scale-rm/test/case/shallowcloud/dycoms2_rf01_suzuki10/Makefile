################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..


# user-defined source files
CODE_DIR    = ..
ORG_SRCS    = mod_user.f90

# parameters for run
INITCONF    = init.conf
RUNCONF     = run.conf
TPROC       = 144

# required data (parameters,distributed files)
DATDIR      = ../../../data
DATPARAM    =
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

analy:
	$(FC) $(FFLAGS_FAST) -o .libs/profile_maker ./postprocess/make_profile_cmp_nc_each.f90 $(NETCDF_INCLUDE) $(NETCDF_LIBS)
