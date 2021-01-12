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
PPNAME      = gm_ll2ico
PPCONF      = pp.conf
TPROC       = 10

# required data (parameters,distributed files)
DATPARAM    = " [$(SCALE_DB)/topo/GTOPO30/Products input_topo] "
DATDISTS    = " [${TPROC} $(SCALE_DB)/grid/boundary/gl05rl00pe10/boundary_GL05RL00 boundary_GL05RL00] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
