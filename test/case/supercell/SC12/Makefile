TOPDIR = ../../../..

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude


index     = SC200m50x50x50
tracer    = kessler
radiation = mstrnX
sfcflx    = louis
precision = real8
postfix   = 

all: build

build:
	cd $(SRCDIR) && $(MAKE) index=$(index) tracer=$(tracer) radiation=$(radiation) sfcflx=$(sfcflx) precision=$(precision) postfix=$(postfix)

runinit: build
	run.sh $(BINDIR)/$(INITNAME) init.conf 1

run: build
	run.sh $(BINDIR)/$(BINNAME) run.conf 1

check:
# to be add program to check result

distclean:
	$(MAKE) -C $(SRCDIR) distclean
