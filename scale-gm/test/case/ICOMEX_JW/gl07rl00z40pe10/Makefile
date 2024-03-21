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
DATPARAM    = " [${DATADIR_GM}/grid/vgrid/vgrid40_stretch_45km.dat vgrid40_stretch_45km.dat] \
                [$(SCALE_DB)/grid/llmap/gl07rl00pe10_i720j360 llmap] "
DATDISTS    = " [${TPROC} $(SCALE_DB)/grid/boundary/gl07rl00pe10/boundary_GL07RL00 boundary_GL07RL00] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
