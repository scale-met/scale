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
TPROC       = 10

# required data (parameters,distributed files)
DATPARAM    = " [${DATADIR_GM}/grid/vgrid/vgrid10_1000m_dcmip.dat vgrid10_1000m_dcmip.dat] \
                [$(SCALE_DB)/grid/llmap/gl05rl00pe10_i180j90 llmap] "
DATDISTS    = " [${TPROC} $(SCALE_DB)/grid/boundary/gl05rl00pe10/boundary_GL05RL00 boundary_GL05RL00] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
