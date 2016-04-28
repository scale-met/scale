################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

# optional parameters for build
index       = 500m200x40x2
precision   = real8
dynamics    = heve
tracer      = dry
sfcflx      = const
turbulence  = dummy
radiation   = dummy

# user-defined source files
CODE_DIR    = ..
ORG_SRCS    = mod_user.f90

# parameters for run
INITCONF    = init.conf
RUNCONF     = run.conf
TPROC       = 150

# required data (parameters,distributed files)
DATDIR      = ../../../data
DATPARAM    =
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
