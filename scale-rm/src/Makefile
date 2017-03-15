################################################################################
#
# Makefile for main program
#
################################################################################

TOPDIR      = $(abspath ../..)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = ../../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

BINNAME  = scale-rm
INITNAME = scale-rm_init
PPNAME   = scale-rm_pp

BINS = $(BINDIR)/$(BINNAME) $(BINDIR)/$(INITNAME) $(BINDIR)/$(PPNAME)
LIBS = $(LIBDIR)/libscale.a $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a

ifneq ($(IGNOREDEP),T)
  DEPENDLIB = $(LIBS)
endif

VPATH = $(BUILD_DIR): \
        admin:        \
        atmos:        \
        ocean:        \
        land:         \
        urban:        \
        coupler:      \
        user:         \
        preprocess:

VERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(VERSION),)
  VERSION  = $(shell cat VERSION)
else
  VERSION := $(VERSION)
endif

SUB_OBJS =	\
	mod_admin_time.o          \
	mod_atmos_admin.o         \
	mod_ocean_admin.o         \
	mod_land_admin.o          \
	mod_urban_admin.o         \
	mod_cpl_admin.o           \
	\
	mod_atmos_dyn_vars.o      \
	mod_atmos_phy_cp_vars.o   \
	mod_atmos_phy_mp_vars.o   \
	mod_atmos_phy_rd_vars.o   \
	mod_atmos_phy_sf_vars.o   \
	mod_atmos_phy_tb_vars.o   \
	mod_atmos_phy_ch_vars.o   \
	mod_atmos_phy_ae_vars.o   \
	mod_atmos_vars.o          \
	mod_ocean_vars.o          \
	mod_land_vars.o           \
	mod_urban_vars.o          \
	mod_cpl_vars.o            \
	\
	mod_admin_restart.o       \
	mod_atmos_dyn_driver.o    \
	mod_atmos_phy_mp_driver.o \
	mod_atmos_phy_ae_driver.o \
	mod_atmos_phy_ch_driver.o \
	mod_atmos_phy_rd_driver.o \
	mod_atmos_phy_sf_driver.o \
	mod_atmos_phy_tb_driver.o \
	mod_atmos_phy_cp_driver.o \
	mod_atmos_driver.o        \
	mod_ocean_phy_driver.o    \
	mod_ocean_driver.o        \
	mod_land_phy_driver.o     \
	mod_land_driver.o         \
	mod_urban_phy_driver.o    \
	mod_urban_driver.o        \
	mod_cpl_driver.o          \
	\
	mod_copytopo.o            \
	mod_cnvtopo.o             \
	mod_cnvlanduse.o          \
	mod_convert.o             \
	mod_mktopo.o              \
	mod_realinput_scale.o     \
	mod_realinput_wrfarw.o    \
	mod_realinput_nicam.o     \
	mod_realinput_grads.o     \
	mod_realinput_scale.o     \
	mod_realinput.o           \
	mod_mkinit.o

MAIN_OBJS =	\
	mod_rm_driver.o           \
	mod_rm_prep.o

OBJS = 	\
	$(SUB_OBJS)	\
	$(MAIN_OBJS)

all:
	$(MAKE) makedir
	$(MAKE) -C $(DCUTILSDIR)
	$(MAKE) -C $(GTOOLDIR)
	$(MAKE) -C $(SCALELIBDIR)/src BUILD_DIR=$(BUILD_DIR)
	@echo;echo "Entering scale-rm...";echo "Current version is " $(VERSION)
	$(MAKE) makebin
	@echo "Complete making scale-rm."

makedir:
	mkdir -p $(BINDIR)
	mkdir -p $(LIBDIR)
	mkdir -p $(BUILD_DIR)

$(LIBDIR)/libdcutils.a :
	$(MAKE) -C $(DCUTILSDIR)

$(LIBDIR)/libgtool.a : $(LIBDIR)/libdcutils.a
	$(MAKE) -C $(GTOOLDIR)

$(LIBDIR)/libscale.a : $(LIBDIR)/libdcutils.a $(LIBDIR)/libgtool.a
	$(MAKE) -C $(SCALELIBDIR)/src

makebin: $(BINS)

$(BINDIR)/$(BINNAME)  : $(BUILD_DIR)/$(BINNAME).o  $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o $(LIBS)
	$(FC) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(NETCDF_LIBS) $(CONTRIB_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)
$(BINDIR)/$(INITNAME) : $(BUILD_DIR)/$(INITNAME).o $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o  $(LIBS)
	$(FC) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(NETCDF_LIBS) $(CONTRIB_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)
$(BINDIR)/$(PPNAME)   : $(BUILD_DIR)/$(PPNAME).o   $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o  $(LIBS)
	$(FC) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(NETCDF_LIBS) $(CONTRIB_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BUILD_DIR)/$(BINNAME).o  : $(BINNAME).f90  $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o
$(BUILD_DIR)/$(INITNAME).o : $(INITNAME).f90 $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o
$(BUILD_DIR)/$(PPNAME).o   : $(PPNAME).f90   $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o

$(BUILD_DIR)/mod_user.o : mod_user.f90 $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(SUB_OBJS))



allclean: distclean
	$(MAKE) -C $(DCUTILSDIR)      allclean
	$(MAKE) -C $(GTOOLDIR)        allclean
	$(MAKE) -C $(SCALELIBDIR)/src allclean
	rm -f $(BINDIR)/$(BINNAME)
	rm -f $(BINDIR)/$(INITNAME)
	rm -f $(BINDIR)/$(PPNAME)

distclean: clean
	rm -f $(BUILD_DIR)/$(BINNAME)
	rm -f $(BUILD_DIR)/$(INITNAME)
	rm -f $(BUILD_DIR)/$(PPNAME)

clean:
	rm -f $(BUILD_DIR)/*.*

dependency:
	./makedepend .



.SUFFIXES:
.SUFFIXES: .o .f90 .c .mod

%.mod: %.f90
	make $(patsubst %.f90,%.o,$<)

$(BUILD_DIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MODDIR) $(CONTRIB_INCLUDE) -I$(SCALELIBDIR)/include -I$(SCALERMDIR)/include $(PAPI_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

$(BUILD_DIR)/%.o : %.c
	$(CC) $(CFLAGS) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

.PHONY : clean distclean allclean dcutils gtool scalelib depend dependency


ifneq ($(IGNOREDEP),T)
-include depend
endif
