################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

# for build
index      = 200x11x11
precision  = real8
dynamics   = heve
tracer     = sn13
sfcflx     = const
turbulence = smg
radiation  = dycoms2

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 6


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
