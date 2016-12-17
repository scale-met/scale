PWD = $(shell pwd)
BUILD_DIR := $(PWD)/.libs

TOPDIR     = ../../..
SYSDEP_DIR = ../../../sysdep

DEBUG=T
export DEBUG
AGGRESSIVE=F
export DEBUG

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

VPATH = $(BUILD_DIR)

LDFLAGS = $(FFLAGS_DEBUG)

PROJ_EXE = unit

#default index & tracer set
index   = 10x10x2
dynamics = heve
tracer  = kessler
radiation = mstrnx
sfcflx  = const
turbulence = smg
precision = real8
postfix = debug

MODS = \
	test_comm.o \
	test_atmos_phy_tb_smg.o \
	test_atmos_dyn.o

all:
	$(MAKE) makedir
	$(MAKE) makelib
	$(MAKE) makebin
	$(MAKE) run

makedir:
	mkdir -p $(BUILD_DIR)

makelib:
	cd $(SCALELIBDIR)/src && $(MAKE) index=$(index) tracer=$(tracer) radiation=$(radiation) sfcflx=$(sfcflx) precision=$(precision) postfix=$(postfix) FFLAGS="$(FFLAGS_DEBUG)" BUILD_DIR=$(BUILD_DIR)

makebin: $(BUILD_DIR)/$(PROJ_EXE)

run:
	mpirun -n 4 $(BUILD_DIR)/$(PROJ_EXE) run.conf

$(BUILD_DIR)/$(PROJ_EXE): $(BUILD_DIR)/$(PROJ_EXE).o $(patsubst %,$(BUILD_DIR)/%,$(MODS)) $(LIBDIR)/libscale.a $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a
	$(LD) $(LDFLAGS) -o $@ $^ $(NETCDF_LIBS) $(LAPACK_LIBS)

$(BUILD_DIR)/test_comm.o: test_comm.f90

$(BUILD_DIR)/test_atmos_phy_tb_smg.o: test_atmos_phy_tb_smg.f90

$(BUILD_DIR)/test_atmos_dyn.o: test_atmos_dyn.f90

$(BUILD_DIR)/unit.o: unit.f90 $(patsubst %,$(BUILD_DIR)/%,$(MODS))

.SUFFIXES: .o .f90
$(BUILD_DIR)/%.o : %.f90
	$(FC) $(FFLAGS_DEBUG) -I$(BUILD_DIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MODDIR) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)

%.o: %.mod

.PHONY: clean distclean

clean:
	rm -f $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.o *~

distclean: clean
	rm -f $(BUILD_DIR)/$(PROJ_EXE)
	rm -f LOG.pe000000
