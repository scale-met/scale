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
PPNAME      = scale-rm_pp
INITNAME    = NONE
BINNAME     = NONE
PPCONF      = pp.conf
INITCONF    = NONE
RUNCONF     = NONE
TPROC       = 16

# required data (parameters,distributed files)
DATDIR      = $(SCALE_DB)/landuse/GLCCv2/Products
DATPARAM    =
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
