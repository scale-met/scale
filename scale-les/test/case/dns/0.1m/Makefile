################################################################################
#
# Makefile for each test program
#
################################################################################

#TOPDIR      = ../../../..
#BUILD_DIR  := $(shell pwd)/.libs
#INDEX_DIR  := $(shell pwd)
#SYSDEP_DIR := $(shell pwd)/$(TOPDIR)/sysdep
PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

#include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
#include $(TOPDIR)/src/Mkinclude

# for build
index      = 500m60x100x100
precision  = real8
dynamics   = heve
tracer     = sn13
tracer     = tomita08
sfcflx     = louis
turbulence = dns
radiation  = dummy

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 16


# src files
CODE_DIR     = ..
ORG_SRCS = mod_user.f90


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

