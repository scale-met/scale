TOPDIR = ../../..

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude


index     = 42x56x56
tracer    = ndw6
radiation = mstrnX
sfcflx    = const
precision = real8
postfix   = 

all: build

build:
	cd $(SRCDIR) && $(MAKE) index=$(index) tracer=$(tracer) radiation=$(radiation) sfcflx=$(sfcflx) precision=$(precision) postfix=$(postfix)

runinit:
	$(SCRDIR)/run.sh $(BINDIR)/$(INITNAME) init.conf 6

run:
	$(SCRDIR)/run.sh $(BINDIR)/$(BINNAME) run.conf 6

check:
# to be add program to check result

distclean:
	$(MAKE) -C $(SRCDIR) distclean
