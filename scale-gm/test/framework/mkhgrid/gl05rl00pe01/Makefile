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
PPNAME      = gm_mkhgrid
PPCONF      = pp.conf
TPROC       = 1

# required data (parameters,distributed files)
DATPARAM    = ""
DATDISTS    = " [${TPROC} ../../mkrawgrid/gl05rl00pe01/rawgrid_GL05RL00 rawgrid_GL05RL00] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
