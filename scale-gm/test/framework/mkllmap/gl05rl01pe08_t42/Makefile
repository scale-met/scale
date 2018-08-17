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
PPNAME      = gm_mkllmap
PPCONF      = pp.conf
ICO2LLNAME  = gm_fio_ico2ll
ICO2LLCONF  = ico2ll.conf
TPROC       = 8

# required data (parameters,distributed files)
DATPARAM    = ""
DATDISTS    = " [${TPROC} ../../mkhgrid/gl05rl01pe08/boundary_GL05RL01 boundary_GL05RL01] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
