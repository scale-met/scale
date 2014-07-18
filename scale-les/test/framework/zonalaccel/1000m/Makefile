################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

# user-defined source files
CODE_DIR    = ..
ORG_SRCS    = mod_user.f90

# for build
index      = 1000m20x25x5
precision  = real8
dynamics   = hevi
tracer     = dry
sfcflx     = const
turbulence = dummy
radiation  = dummy

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 4


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
