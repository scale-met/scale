################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

# for build
index      = 500m40x20x4
precision  = real8
dynamics   = heve
tracer     = sn14
sfcflx     = cpl
turbulence = smg
radiation  = mstrnx

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 8


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
