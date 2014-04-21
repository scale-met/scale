################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..
SYSDEP_DIR  = $(abspath ..)

# for build
index      = 005m1256x32x32
precision  = real8
dynamics   = heve
tracer     = sn14
sfcflx     = louis
turbulence = smg
radiation  = mstrnx

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 6


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
