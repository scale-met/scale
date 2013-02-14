TOPDIR = ../../../..
BUILD_DIR := $(shell pwd)/.libs
MKPARADIR = ../../../../src/preprocess/mk_para2
MKPARAFILE = micpara.dat

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

VPATH = $(BUILD_DIR)

# for build
index     = 500m40x40x40
tracer    = hbinw
radiation = mstrnX
sfcflx    = const
precision = real8
postfix   = 

# for run
INITCONF  = init.conf
RUNCONF   = run.conf

all: makedir build buildmkpara

makedir:
	mkdir -p $(BUILD_DIR)

build: makedir
	$(MAKE) -C $(SRCDIR) index=$(index) tracer=$(tracer) radiation=$(radiation) sfcflx=$(sfcflx) precision=$(precision) postfix=$(postfix) BUILD_DIR=$(BUILD_DIR) INDEX_DIR=$(shell pwd)

buildmkpara:
	rm -rf $(MKPARADIR)/mkpara $(MKPARADIR)/*.o $(MKPARADIR)/*.mod
	cp ./setup.h $(MKPARADIR)/
	make -C $(MKPARADIR)
	cp $(MKPARADIR)/mkpara $(BUILD_DIR)
	cp $(MKPARADIR)/*asc* ./ 
	$(BUILD_DIR)/mkpara
	rm -rf ./*asc*

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
	rm -f $(MKPARADIR)/mkpara $(MKPARADIR)/*.o $(MKPARADIR)/*.mod micpara.dat
distclean:
	rm -f *.nc LOG.pe*
