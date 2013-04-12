TOPDIR = ../../../..
BUILD_DIR := $(shell pwd)/.libs

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

VPATH = $(BUILD_DIR)

# for build
index     = 100m80x15x15
tracer    = sn13
radiation = mstrnX
sfcflx    = const
precision = real8
turbulence = dummy
postfix   = 

# for run
INITCONF  = init.conf
RUNCONF   = run.conf

all: build

$(BUILD_DIR)/mod_user.f90: mod_user.f90
	mkdir -p $(BUILD_DIR)
	cp -f $< $@

build: $(BUILD_DIR)/mod_user.f90
	$(MAKE) -C $(SRCDIR) index=$(index) tracer=$(tracer) radiation=$(radiation) sfcflx=$(sfcflx) precision=$(precision) turbulence=$(turbulence) postfix=$(postfix) BUILD_DIR=$(BUILD_DIR) INDEX_DIR=$(shell pwd)

jobshell: build
	sh $(SCRDIR)/Mkjobshell.sh $(BINDIR) ./$(INITNAME) ./$(BINNAME) $(INITCONF) $(RUNCONF) 8

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
