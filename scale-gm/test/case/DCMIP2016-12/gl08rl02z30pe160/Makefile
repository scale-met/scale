################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..


# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
RUNCONF     = run.conf
TPROC       = 160

# required data (parameters,distributed files)
DATPARAM    = " [${TESTDIR}/data/grid/vgrid/vgrid30_stretch_30km_dcmip2016v2.dat vgrid30_stretch_30km_dcmip2016v2.dat] "
DATDISTS    = " [${TPROC} ${TESTDIR}/data/grid/boundary/gl08rl02pe160/boundary_GL08RL02 boundary_GL08RL02] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
