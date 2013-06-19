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
index      = 500m60x100x100
precision  = real8
dyanamics  = heve
tracer     = tomita08
sfcflx     = const
turbulence = dummy
radiation  = dummy
postfix    = 

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 6



all: buildlocal build

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common



jobshell:
	sh ../Mkjobshell.$(SCALE_SYS).sh $(BINDIR) $(INITNAME) $(BINNAME) $(INITCONF) $(RUNCONF) $(TPROC)

buildlocal: makedir $(BUILD_DIR)/mod_user.f90

$(BUILD_DIR)/mod_user.f90: ../mod_user.f90
	cp -f $< $@
