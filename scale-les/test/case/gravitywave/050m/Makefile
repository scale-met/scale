################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

# for build
index      = 50m200x40x2
precision  = real8
dynamics   = heve
tracer     = dry
sfcflx     = const
turbulence = dummy
radiation  = dummy

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = run.conf


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
