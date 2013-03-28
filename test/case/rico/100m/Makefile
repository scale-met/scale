TOPDIR = ../../../..
BUILD_DIR := $(shell pwd)/.libs
CODE_DIR := $(shell pwd)/code

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

VPATH = $(BUILD_DIR)

# for build
index     = 100x16x16
tracer    = ndw6
radiation = rico
sfcflx    = rico
precision = real8
postfix   = 

# for run
INITCONF  = init.conf
RUNCONF   = run.conf

all: makedir build

makedir:

$(BUILD_DIR)/mod_atmos_phy_sf_rico.f90: $(CODE_DIR)/mod_atmos_phy_sf_rico.f90
	mkdir -p $(BUILD_DIR)
	cp -f $< $@
$(BUILD_DIR)/mod_atmos_phy_rd_rico.f90: $(CODE_DIR)/mod_atmos_phy_rd_rico.f90
	mkdir -p $(BUILD_DIR)
	cp -f $< $@
$(BUILD_DIR)/mod_mkinit.f90: $(CODE_DIR)/mod_mkinit_rico.f90
	mkdir -p $(BUILD_DIR)
	cp -f $< $@

build: $(BUILD_DIR)/mod_atmos_phy_sf_rico.f90 $(BUILD_DIR)/mod_atmos_phy_rd_rico.f90 $(BUILD_DIR)/mod_mkinit.f90 
	$(MAKE) -C $(SRCDIR) index=$(index) tracer=$(tracer) radiation=$(radiation) sfcflx=$(sfcflx) precision=$(precision) postfix=$(postfix) BUILD_DIR=$(BUILD_DIR) INDEX_DIR=$(shell pwd)

jobshell: build
	sh $(SCRDIR)/Mkjobshell.sh $(BINDIR) $(INITNAME) $(BINNAME) $(INITCONF) $(RUNCONF) 6

run: jobshell
	$(JOBSUB) run.sh

check:
# to be add program to check result

.PHONY: clean distclean

#allclean:
#	$(MAKE) -C $(SRCDIR) allclean
clean:
	rm -f $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.o *~
distclean:
	rm -f *.nc LOG.pe*
