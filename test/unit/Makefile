TOPDIR = ../..
BUILD_DIR := $(shell pwd)/.libs

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

VPATH = $(BUILD_DIR)

LDFLAGS = $(FFLAGS_DEBUG)

PROJ_EXE = unit

#default index & tracer set
index   = 10x10x2
dynamics = heve
tracer  = kessler
radiation = mstrnX
sfcflx  = const
turbulence = smg
precision = real8
postfix = debug

MODS = \
	test_atmos_phy_tb_smg.o \
	test_atmos_dyn_fent_fct.o

all: makedir makelib makeinclude makebin run

makedir:
	mkdir -p $(BUILD_DIR)

makelib:
	cd $(SRCDIR) && $(MAKE) INDEX_DIR=$(shell pwd) index=$(index) tracer=$(tracer) radiation=$(radiation) sfcflx=$(sfcflx) precision=$(precision) postfix=$(postfix) FFLAGS="$(FFLAGS_DEBUG)" BUILD_DIR=$(BUILD_DIR)

makeinclude:
	cp $(INCDIR)/inc_tracer_$(tracer).f90       $(BUILD_DIR)/inc_tracer.h
	cp $(INCDIR)/inc_precision_$(precision).f90 $(BUILD_DIR)/inc_precision.h

makebin: $(BUILD_DIR)/$(PROJ_EXE)

run:
	$(BUILD_DIR)/$(PROJ_EXE) dummy.conf

$(BUILD_DIR)/$(PROJ_EXE): $(BUILD_DIR)/$(PROJ_EXE).o $(patsubst %,$(BUILD_DIR)/%,$(MODS)) $(LIBDIR)/$(LIBNAME) $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a
	$(LD) $(LDFLAGS) -o $@ $^ $(NETCDF_LIBS)

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
