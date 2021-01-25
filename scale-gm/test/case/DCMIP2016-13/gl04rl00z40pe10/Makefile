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
DATPARAM    = " [${TESTDIR}/data/grid/vgrid/vgrid40_cnst_20km_dcmip2016.dat vgrid40_cnst_20km_dcmip2016.dat] \
                [$(SCALE_DB)/grid/llmap/gl04rl00pe10_i90j45 llmap] "
DATDISTS    = " [${TPROC} $(SCALE_DB)/grid/boundary/gl04rl00pe10/boundary_GL04RL00 boundary_GL04RL00] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
