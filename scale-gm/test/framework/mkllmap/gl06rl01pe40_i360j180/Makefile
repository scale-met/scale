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
TPROC       = 40

# required data (parameters,distributed files)
DATPARAM    = ""
DATDISTS    = " [${TPROC} ../../mkhgrid/gl06rl01pe40/boundary_GL06RL01 boundary_GL06RL01] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
