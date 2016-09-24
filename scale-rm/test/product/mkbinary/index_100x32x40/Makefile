################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..
INDEX_DIR   = .

# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# for build
index      = 100x32x40
precision  = real8
dynamics   = hevi
tracer     = tomita08
sfcflx     = couple
turbulence = smagorinsky
radiation  = mstrnx

# parameters for run
INITCONF    =
RUNCONF     =
TPROC       =

# required data (parameters,distributed files)
DATDIR      =
DATPARAM    =
DATDISTS    =

# options
export AGGRESSIVE        = T
export ENABLE_FIXEDINDEX = T

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
