TOPDIR = ../../..
BUILD_DIR := $(shell pwd)/.libs

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

VPATH = $(BUILD_DIR)

index     = 30m10a82x16x16
tracer    = dry
radiation = mstrnX
sfcflx    = const
precision = real8
postfix   = 

all: build

build:
	cd $(SRCDIR) && $(MAKE) index=$(index) tracer=$(tracer) radiation=$(radiation) sfcflx=$(sfcflx) precision=$(precision) postfix=$(postfix) BUILD_DIR=$(BUILD_DIR) INDEX_DIR=$(shell pwd)

runinit: build
	$(SCRDIR)/run.sh $(BINDIR)/$(INITNAME) init_30m10a.conf 4

run: build
	$(SCRDIR)/run.sh $(BINDIR)/$(BINNAME) run_30m10a.conf 4

check:
# to be add program to check result

.PHONY: clean distclean

clean:
	rm -f $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.o *~
distclean:
	rm -f *.nc LOG.pe000000
