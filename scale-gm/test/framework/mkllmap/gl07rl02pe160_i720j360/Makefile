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
PPNAME      = gm_mkllmap$(POSTFIX)
PPCONF      = pp.conf
ICO2LLCONF  = ico2ll.conf
TPROC       = 160

# required data (parameters,distributed files)
DATPARAM    = ""
DATDISTS    = " [${TPROC} ../../mkhgrid/gl07rl02pe160/boundary_GL07RL02 boundary_GL07RL02] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
