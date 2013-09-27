################################################################################
#
# Makefile for each test program
#
################################################################################

TOPDIR      = ../../../..
BUILD_DIR  := $(shell pwd)/.libs
INDEX_DIR  := $(shell pwd)
SYSDEP_DIR := $(shell pwd)/$(TOPDIR)/sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

# for build
index      = 500m40x40x40
dynamics   = heve
precision  = real8
tracer     = sn13
sfcflx     = const
turbulence = dummy
radiation  = mstrnX
postfix    = 

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 6



all: build

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

