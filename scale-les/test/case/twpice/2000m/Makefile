################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

# for build
index      = 2000m120x25x25
precision  = real8
dynamics   = heve
tracer     = tomita08
sfcflx     = louis
turbulence = smg
radiation  = dummy

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 16


# src files
CODE_DIR = ..
ORG_SRCS = mod_user.f90


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

