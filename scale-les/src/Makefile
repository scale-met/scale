################################################################################
#
# Makefile for main program
#
################################################################################

TOPDIR      = $(abspath ../..)
BUILD_DIR   = ./.libs
INDEX_DIR   = ./include
SYSDEP_DIR  = ../sysdep

include $(TOPDIR)/Mkinclude
include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)

MAKEOVERRIDES =

# for build
index        = 500m40x40x40
precision    = real8
dynamics     = heve
tracer       = kessler
sfcflx       = louis
turbulence   = smg
microphysics = $(tracer)
radiation    = dummy
aerosol      = dummy
chemistry    = dummy
postfix      =

VPATH = $(BUILD_DIR):include:atmos:land:ocean:cpl:preprocess:user

VERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(VERSION),)
  VERSION = $(shell cat VERSION)
else
  VERSION := $(VERSION)
endif



PRJ1 = scale3
PRJ2 = scale3_init

PRJS = $(PWD)/$(PRJ1) $(PWD)/$(PRJ2)

MODS =	\
	mod_atmos_vars.o	\
	mod_atmos_vars_sf.o	\
	mod_land_vars.o	\
	mod_ocean_vars.o	\
	mod_cpl_vars.o \
	mod_atmos_dyn_driver.o	\
	mod_atmos_phy_mp_driver.o	\
	mod_atmos_phy_sf_driver.o	\
	mod_atmos_phy_tb_driver.o	\
	mod_atmos_phy_ae_driver.o \
	mod_atmos_phy_rd_driver.o \
	mod_cpl_atmos_land.o \
	mod_cpl.o \
	mod_atmos.o	\
	mod_land_phy_bucket.o	\
	mod_land.o	\
	mod_ocean_phy_fixedsst.o	\
	mod_ocean.o	\
	mod_mktopo.o	\
	mod_mkinit.o	\
	mod_user.o

all: makedir makedcutils makegtool makescalelib makeinclude makebin

makedcutils:
	cd $(DCUTILSDIR); $(MAKE)

makegtool:
	cd $(GTOOLDIR); $(MAKE)

makescalelib:
	cd $(SCALELIBDIR)/src; $(MAKE) depend
	cd $(SCALELIBDIR)/src; $(MAKE)

makedir:
	mkdir -p $(BINDIR)
	mkdir -p $(BUILD_DIR)

makeinclude:
	cp $(INDEX_DIR)/inc_land.f90 $(BUILD_DIR)/inc_land.h

makebin: $(PRJS)

$(PWD)/$(PRJ1): $(BUILD_DIR)/$(PRJ1).o $(patsubst %,$(BUILD_DIR)/%,$(MODS))
	$(FC) $(LDFLAGS) -o $@ $^  $(LIBDIR)/libscalelib.a $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(PWD)/$(PRJ2): $(BUILD_DIR)/$(PRJ2).o $(patsubst %,$(BUILD_DIR)/%,$(MODS))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBDIR)/libscalelib.a $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BUILD_DIR)/$(PRJ1).o: $(PRJ1).f90 $(patsubst %,$(BUILD_DIR)/%,$(MODS))
$(BUILD_DIR)/$(PRJ2).o: $(PRJ2).f90 $(patsubst %,$(BUILD_DIR)/%,$(MODS))

$(BUILD_DIR)/mod_atmos.o	: mod_atmos.f90 mod_atmos_vars.o \
                             mod_atmos_dyn_driver.o mod_atmos_phy_sf_driver.o mod_atmos_phy_tb_driver.o mod_atmos_phy_rd_driver.o mod_atmos_phy_ae_driver.o mod_cpl_vars.o
$(BUILD_DIR)/mod_atmos_vars.o	: mod_atmos_vars.f90
$(BUILD_DIR)/mod_atmos_vars_sf.o	: mod_atmos_vars_sf.f90
$(BUILD_DIR)/mod_atmos_dyn_driver.o	: mod_atmos_dyn_driver.f90 mod_atmos_vars.o mod_atmos_vars_sf.o mod_ocean_vars.o
$(BUILD_DIR)/mod_atmos_phy_sf_driver.o	: mod_atmos_phy_sf_driver.f90 mod_atmos_vars.o mod_atmos_vars_sf.o mod_ocean_vars.o
$(BUILD_DIR)/mod_atmos_phy_tb_driver.o	: mod_atmos_phy_tb_driver.f90 mod_atmos_vars.o
$(BUILD_DIR)/mod_atmos_phy_rd_driver.o	: mod_atmos_phy_rd_driver.f90 mod_atmos_phy_ae_driver.o mod_atmos_vars.o mod_atmos_vars_sf.o mod_cpl_vars.o

$(BUILD_DIR)/mod_land.o	: mod_land.f90 mod_land_vars.o mod_land_phy_bucket.o mod_cpl_vars.o
$(BUILD_DIR)/mod_land_vars.o	: mod_land_vars.f90
$(BUILD_DIR)/mod_land_phy_bucket.o	: mod_land_phy_bucket.f90 mod_land_vars.o

$(BUILD_DIR)/mod_ocean.o	: mod_ocean.f90 mod_ocean_vars.o mod_ocean_phy_fixedsst.o
$(BUILD_DIR)/mod_ocean_vars.o	: mod_ocean_vars.f90
$(BUILD_DIR)/mod_ocean_phy_fixedsst.o	: mod_ocean_phy_fixedsst.f90 mod_ocean_vars.o

$(BUILD_DIR)/mod_cpl.o	: mod_cpl.f90 mod_cpl_vars.o mod_cpl_atmos_land.o
$(BUILD_DIR)/mod_cpl_vars.o	: mod_cpl_vars.f90
$(BUILD_DIR)/mod_cpl_atmos_land.o	: mod_cpl_atmos_land.f90 mod_cpl_vars.o

$(BUILD_DIR)/mod_mktopo.o	: mod_mktopo.f90
$(BUILD_DIR)/mod_mkinit.o	: mod_mkinit.f90 mod_atmos_vars.o

$(BUILD_DIR)/mod_user.o	: mod_user.f90

$(BUILD_DIR)/scale3_init.o	: scale3_init.f90 mod_atmos_vars.o mod_mktopo.o mod_mkinit.o
$(BUILD_DIR)/scale.o			: scale.f90 mod_atmos.o mod_atmos_vars.o mod_ocean.o mod_ocean_vars.o mod_cpl.o mod_cpl_vars.o

allclean: distclean
	$(MAKE) -C $(DCUTILSDIR) allclean
	$(MAKE) -C $(GTOOLDIR) allclean
	$(MAKE) -C $(SCALELIBDIR)/src allclean
clean:
	rm -f $(PRJS) $(BUILD_DIR)/*.h $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.o $(BUILD_DIR)/*.lst

distclean: clean
	rm -f $(PRJS)

.SUFFIXES:
.SUFFIXES: .o .f90 .c

$(BUILD_DIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INDEX_DIR) -I$(SCALELIBDIR)/include -I$(INCDIR) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/mod_atmos_dyn.o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INDEX_DIR) -I$(INCDIR) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/mod_atmos_dyn_$(dynamics).o:
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INDEX_DIR) -I$(INCDIR) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/mod_atmos_phy_rd_profile.o:
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(GTOOLDIR) -I$(DCUTILSDIR) -I$(INDEX_DIR) -I$(INCDIR) $(NETCDF_INCLUDE) $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)


.PHONY : clean distclean allclean
