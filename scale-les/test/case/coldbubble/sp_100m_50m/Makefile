################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

# for build
index      = 100m64x128x2
precision  = real8
dynamics   = heve
tracer     = dry
sfcflx     = const
turbulence = dummy
radiation  = dummy

# for run
INITCONF  = init_100m.conf
RUNCONF   = run_100m.conf
TPROC     = 4


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
