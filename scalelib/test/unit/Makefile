PWD = $(shell pwd)
TOPDIR = ../../..
BUILD_DIR := $(PWD)/.libs

include $(TOPDIR)/Mkinclude
include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)

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
	test_atmos_phy_tb_smg.o \
	test_atmos_dyn_fent_fct.o

all: makedir makelib makebin run

makedir:
	mkdir -p $(BUILD_DIR)

makelib:
	cd $(SCALELIBDIR)/src && $(MAKE) index=$(index) tracer=$(tracer) radiation=$(radiation) sfcflx=$(sfcflx) precision=$(precision) postfix=$(postfix) FFLAGS="$(FFLAGS_DEBUG)" BUILD_DIR=$(BUILD_DIR)

makebin: $(BUILD_DIR)/$(PROJ_EXE)

run:
	$(BUILD_DIR)/$(PROJ_EXE) dummy.conf

$(BUILD_DIR)/$(PROJ_EXE): $(BUILD_DIR)/$(PROJ_EXE).o $(patsubst %,$(BUILD_DIR)/%,$(MODS)) $(LIBDIR)/libscale.a $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a
	$(LD) $(LDFLAGS) -o $@ $^ $(NETCDF_LIBS) $(LAPACK_LIBS)

$(BUILD_DIR)/test_atmos_phy_tb_smg.o: test_atmos_phy_tb_smg.f90

$(BUILD_DIR)/test_atmos_dyn_fent_fct.o: test_atmos_dyn_fent_fct.f90

$(BUILD_DIR)/unit.o: unit.f90 test_atmos_phy_tb_smg.o test_atmos_dyn_fent_fct.o

.SUFFIXES: .o .f90
$(BUILD_DIR)/%.o : %.f90
	$(FC) $(FFLAGS_DEBUG) -I$(BUILD_DIR) -I$(DCUTILSDIR) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)

%.o: %.mod

.PHONY: clean distclean

clean:
	rm -f $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.o *~

distclean: clean
	rm -f $(BUILD_DIR)/$(PROJ_EXE)
	rm -f LOG.pe000000
