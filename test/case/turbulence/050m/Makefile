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
index      = 050m120x32x32
precision  = real8
dynamics   = heve
tracer     = dry
sfcflx     = const
turbulence = smg
radiation  = dummy
postfix    = 

ADDITIONAL_FFLAGS = -DDRY

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 150



all: build

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
