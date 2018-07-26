################################################################################
#
# Makefile for library unit test
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../..)
SYSDEP_DIR  = $(TOPDIR)/sysdep
BUILD_DIR  := $(PWD)/.libs
BINDIR      = $(PWD)

SCALE_DEBUG=T
export SCALE_DEBUG
SCALE_USE_AGRESSIVEOPT=F
export SCALE_USE_AGRESSIVEOPT

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

VPATH = $(BUILD_DIR)

BINNAME = unit

LIBS = $(LIBDIR)/libscale.a $(LIBDIR)/libdcutils.a

OBJS = \
	test_comm.o             \
	test_atmos_phy_tb_smg.o \
	test_atmos_dyn.o

all:
	$(MAKE) makedir
	$(MAKE) makelib
	$(MAKE) makebin
	$(MAKE) run

makedir:
	mkdir -p $(LIBDIR)
	mkdir -p $(BUILD_DIR)

makelib:
	$(MAKE) -C $(DCUTILSDIR)
	$(MAKE) -C $(SCALELIBDIR)/src BUILD_DIR=$(BUILD_DIR)

makebin: $(BINDIR)/$(BINNAME)

run:
	mpirun -n 4 ./$(BINNAME) run.conf

$(BINDIR)/$(BINNAME)  : $(BUILD_DIR)/$(BINNAME).o  $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS) $(CONTRIB_LIBS)

$(BUILD_DIR)/$(BINNAME).o  : $(BINNAME).f90  $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(OBJS))

$(BUILD_DIR)/test_comm.o:             test_comm.f90
$(BUILD_DIR)/test_atmos_phy_tb_smg.o: test_atmos_phy_tb_smg.f90
$(BUILD_DIR)/test_atmos_dyn.o:        test_atmos_dyn.f90
$(BUILD_DIR)/unit.o:                  unit.f90 $(patsubst %,$(BUILD_DIR)/%,$(MODS))



allclean: distclean
	$(MAKE) -C $(DCUTILSDIR)      allclean
	$(MAKE) -C $(SCALELIBDIR)/src allclean
	rm -f $(BINDIR)/$(BINNAME)

distclean: clean
	rm -f $(BUILD_DIR)/$(BINNAME)
	rm -f LOG.pe000000

clean:
	rm -rf $(BUILD_DIR)/*.dSYM
	rm -f  $(BUILD_DIR)/*.*



.SUFFIXES:
.SUFFIXES: .o .f90 .mod

%.mod: %.f90
	make $(patsubst %.f90,%.o,$<)

$(BUILD_DIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(DCUTILSDIR) -I$(MODDIR) -I$(SCALELIBDIR)/include -I$(SCALERMDIR)/include $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

.PHONY : clean distclean allclean dcutils scalelib

include $(TOPDIR)/utils/make/Make_environments
