################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

# for build
index      = 1000m40x60x60
precision  = real8
dynamics   = heve
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
