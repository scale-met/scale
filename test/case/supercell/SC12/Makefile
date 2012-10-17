TOPDIR = ../../../..
BUILD_DIR := $(shell pwd)/.libs

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

VPATH = $(BUILD_DIR)

index     = SC200m50x50x50
tracer    = kessler
radiation = mstrnX
sfcflx    = louis
precision = real8
postfix   = 

all: build

build:
	cd $(SRCDIR) && $(MAKE) index=$(index) tracer=$(tracer) radiation=$(radiation) sfcflx=$(sfcflx) precision=$(precision) postfix=$(postfix) BUILD_DIR=$(BUILD_DIR)

runinit: build
	$(SCRDIR)/run.sh $(BINDIR)/$(INITNAME) init.conf 4

run: build
	$(SCRDIR)/run.sh $(BINDIR)/$(BINNAME) run.conf 4

check:
# to be add program to check result

.PHONY: clean

clean:
	rm -f $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.o *~
distclean:
	rm -f *.nc LOG.pe000000
