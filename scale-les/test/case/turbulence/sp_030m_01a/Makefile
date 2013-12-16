################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

# for build
index      = 30m1a82x16x16
precision  = real8
dynamics   = heve
tracer     = dry
sfcflx     = const
turbulence = smg
radiation  = dummy

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = run.conf


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
