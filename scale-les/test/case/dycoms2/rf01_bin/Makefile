################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

# for build
index      = 200x8x8
precision  = real8
dynamics   = heve
tracer     = binw
sfcflx     = const
turbulence = smg
radiation  = dycoms2

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = run.conf


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
