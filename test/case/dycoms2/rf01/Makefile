TOPDIR = ../../../..
BUILD_DIR := $(shell pwd)/.libs
CODE_DIR := $(shell pwd)/code

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

VPATH = $(BUILD_DIR)

# for build
index     = 200x8x8
tracer    = ndw6
radiation = dycoms2
sfcflx    = const
precision = real8
postfix   = 

# for run
INITCONF  = init.conf
RUNCONF   = run.conf

all: makedir build analy

makedir:

$(BUILD_DIR)/mod_atmos_phy_rd_dycoms2.f90: $(CODE_DIR)/mod_atmos_phy_rd_dycoms2.f90
	mkdir -p $(BUILD_DIR)
	cp -f $< $@

build: $(BUILD_DIR)/mod_atmos_phy_rd_dycoms2.f90 
	$(MAKE) -C $(SRCDIR) index=$(index) tracer=$(tracer) radiation=$(radiation) sfcflx=$(sfcflx) precision=$(precision) postfix=$(postfix) BUILD_DIR=$(BUILD_DIR) INDEX_DIR=$(shell pwd)

jobshell: build
	sh $(SCRDIR)/Mkjobshell.sh $(BINDIR) $(INITNAME) $(BINNAME) $(INITCONF) $(RUNCONF) 6

run: jobshell
	$(JOBSUB) run.sh

analy: 
	$(FC) $(FFLAGS) -L$(NETCDF_LIBS) -lnetcdff -lnetcdf -I$(NETCDF_INCLUDE) -o $(BUILD_DIR)/profile_maker $(CODE_DIR)/make_profile_cmp_nc_each.f90
check:
# to be add program to check result

.PHONY: clean distclean

allclean:
	$(MAKE) -C $(SRCDIR) allclean
clean:
	rm -f $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.o *~ $(BUILD_DIR)/profile_maker
distclean:
	rm -f *.nc LOG.pe*
