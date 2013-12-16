################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

# for build
index      = 500m40x40x40
precision  = real8
dynamics   = heve
tracer     = sn13
sfcflx     = const
turbulence = dummy
radiation  = mstrnX

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = run.conf


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
