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
ICO2LLNAME  = gm_fio_ico2ll
ICO2LLCONF  = ico2ll.conf
TPROC       = 10

# required data (parameters,distributed files)
DATPARAM    = " [${TESTDIR}/data/grid/vgrid/vgrid60_200m_dcmip.dat vgrid60_200m_dcmip.dat] \
                [$(SCALE_DB)/grid/llmap/gl06rl00pe10_i360j180 llmap] "
DATDISTS    = " [${TPROC} $(SCALE_DB)/grid/boundary/gl06rl00pe10/boundary_GL06RL00 boundary_GL06RL00] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
