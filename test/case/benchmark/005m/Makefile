TOPDIR = ../../../..
BUILD_DIR := $(shell pwd)/.libs

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

VPATH = $(BUILD_DIR)

# for build
index      = 005m1256x32x32
precision  = real8
tracer     = ndw6
sfcflx     = const
radiation  = dummy
turbulence = dummy
postfix    = 

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 2x3

all: makedir build

makedir:
	mkdir -p $(BUILD_DIR)

build: makedir
	$(MAKE) -C $(SRCDIR) INDEX_DIR=$(shell pwd)   \
                        BUILD_DIR=$(BUILD_DIR)   \
                        index=$(index)           \
                        precision=$(precision)   \
                        tracer=$(tracer)         \
                        sfcflx=$(sfcflx)         \
                        radiation=$(radiation)   \
                        turbulence=$(turbulence) \
                        postfix=$(postfix)

jobshell:
	sh $(TESTDIR)/sysdep/Mkjobshell.$(SCALE_SYS).sh $(BINDIR) $(INITNAME) $(BINNAME) $(INITCONF) $(RUNCONF) $(TPROC)

run: jobshell
	$(JOBSUB) run.sh

check:
# to be add program to check result

.PHONY: allclean clean distclean

allclean: clean
	$(MAKE) -C $(SRCDIR) allclean

clean:
	rm -f $(BUILD_DIR)/*.h $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.o $(BUILD_DIR)/*.a

distclean:
	rm -f *.nc LOG.pe*
