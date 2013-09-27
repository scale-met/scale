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
index      = 400m16x32x2
precision  = real8
dynamics   = heve
tracer     = dry
sfcflx     = const
turbulence = dummy
radiation  = dummy
postfix    = 

ADDITIONAL_FFLAGS = -DDRY

# for run
INITCONF  = init_400m.conf
RUNCONF   = run_400m.conf
TPROC     = 4



all: build

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
