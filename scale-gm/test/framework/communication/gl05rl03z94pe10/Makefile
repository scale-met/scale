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
TPROC       = 10

# required data (parameters,distributed files)
DATPARAM    = " [${TESTDIR}/data/grid/vgrid/vgrid94.dat vgrid94.dat] "
DATDISTS    = " [${TPROC} $(SCALE_DB)/grid/boundary/gl05rl03pe10/boundary_GL05RL03 boundary_GL05RL03] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
