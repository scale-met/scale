################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..
DATADIR     = ${TOPDIR}/data
DATADIR_GM  = ${TESTDIR}/data


# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
RUNCONF     = run.conf
ICO2LLCONF  = ico2ll.conf
TPROC       = 40

# required data (parameters,distributed files)
DATPARAM    = " [${DATADIR_GM}/grid/vgrid/vgrid30_stretch_30km_dcmip2016v2.dat vgrid30_stretch_30km_dcmip2016v2.dat] \
                [$(SCALE_DB)/grid/llmap/gl05rl01pe40_i180j90 llmap] "
DATDISTS    = " [${TPROC} $(SCALE_DB)/grid/boundary/gl05rl01pe40/boundary_GL05RL01 boundary_GL05RL01] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
