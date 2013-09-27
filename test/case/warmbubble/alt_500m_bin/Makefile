################################################################################
#
# Makefile for each test program
#
################################################################################

TOPDIR      = ../../../..
BUILD_DIR  := $(shell pwd)/.libs
MKPARADIR = ../../../../src/preprocess/mk_para2
MKPARAFILE = micpara.dat
INDEX_DIR  := $(shell pwd)
SYSDEP_DIR := $(shell pwd)/$(TOPDIR)/sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

# for build
index      = 500m40x40x40
precision  = real8
dynamics   = heve
tracer     = hbinw
sfcflx     = const
turbulence = smg
radiation  = dummy
postfix    = 

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 6



all: build buildmkpara

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common



buildmkpara:
	rm -rf $(MKPARADIR)/mkpara $(MKPARADIR)/*.o $(MKPARADIR)/*.mod
	cp ./setup.h $(MKPARADIR)/
	make -C $(MKPARADIR)
	cp $(MKPARADIR)/mkpara $(BUILD_DIR)
	cp $(MKPARADIR)/*asc* ./ 
	$(BUILD_DIR)/mkpara
	rm -rf ./*asc*
