TOPDIR = ../../../..
BUILD_DIR := $(shell pwd)/.libs

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

VPATH = $(BUILD_DIR)

# for build
index     = 500m40x40x40
tracer    = kessler
radiation = mstrnX
sfcflx    = const
precision = real8
postfix   = 

# for run
INITCONF  = init.conf
RUNCONF   = run.conf

all: makedir build

makedir:
	mkdir -p $(BUILD_DIR)

build: makedir
	$(MAKE) -C $(SRCDIR) index=$(index) tracer=$(tracer) radiation=$(radiation) sfcflx=$(sfcflx) precision=$(precision) postfix=$(postfix) BUILD_DIR=$(BUILD_DIR) INDEX_DIR=$(shell pwd)

jobshell: build
	sh $(SCRDIR)/Mkjobshell.sh $(BINDIR) $(INITNAME) $(BINNAME) $(INITCONF) $(RUNCONF) 6

run: jobshell
	$(JOBSUB) run.sh

check:
# to be add program to check result

.PHONY: clean distclean

allclean:
	$(MAKE) -C $(SRCDIR) allclean
clean:
	rm -f $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.o *~
distclean:
	rm -f *.nc LOG.pe*
