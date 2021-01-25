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
ICO2LLCONF  = ico2ll.conf
TPROC       = 10

# required data (parameters,distributed files)
DATPARAM    = " [${TESTDIR}/data/grid/vgrid/vgrid120_100m_dcmip.dat vgrid120_100m_dcmip.dat] \
                [$(SCALE_DB)/grid/llmap/gl07rl00pe10_i720j360 llmap] "
DATDISTS    = " [${TPROC} $(SCALE_DB)/grid/boundary/gl07rl00pe10/boundary_GL07RL00 boundary_GL07RL00] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
