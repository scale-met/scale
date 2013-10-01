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
index      = 500m48x8x8
precision  = real8
dynamics   = heve
tracer     = dry
sfcflx     = cpl
turbulence = smg
radiation  = dummy
postfix    = 

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 1



all: build

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
