TOPDIR = ../../..
BUILD_DIR := $(shell pwd)/.libs

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

VPATH = $(BUILD_DIR)

index     = 500m10x56x56
tracer    = ndw6
radiation = mstrnX
sfcflx    = const
precision = real8
postfix   = 

all: build

build:
	cd $(SRCDIR) && $(MAKE) index=$(index) tracer=$(tracer) radiation=$(radiation) sfcflx=$(sfcflx) precision=$(precision) postfix=$(postfix) BUILD_DIR=$(BUILD_DIR) INDEX_DIR=$(shell pwd)

runinit: build
	$(SCRDIR)/run.sh $(BINDIR)/$(INITNAME) init.conf 6

run: build
	$(SCRDIR)/run.sh $(BINDIR)/$(BINNAME) run.conf 6

check:
# to be add program to check result

.PHONY: clean distclean

clean:
	rm -f $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.o *~
distclean:
	rm -f *.nc LOG.pe000000
