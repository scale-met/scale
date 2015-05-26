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

BINNAME  = scale-les
INITNAME = scale-les_init
PPNAME   = scale-les_pp

ifeq ($(NOINIT),)
  BINS = $(BINDIR)/$(BINNAME) $(BINDIR)/$(INITNAME) $(BINDIR)/$(PPNAME)
else
  BINS = $(BINDIR)/$(BINNAME)
endif

DEPENDLIB = $(LIBDIR)/libscale.a $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a

VPATH = $(BUILD_DIR): \
        include:      \
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

OBJS =	\
	mod_admin_restart.o          \
	mod_admin_time.o          \
	mod_atmos_admin.o            \
	mod_land_admin.o             \
	mod_ocean_admin.o            \
	mod_urban_admin.o            \
	mod_cpl_admin.o              \
	mod_cpl_vars.o               \
	mod_cpl_driver.o             \
	mod_atmos_dyn_vars.o         \
	mod_atmos_phy_mp_vars.o      \
	mod_atmos_phy_ae_vars.o      \
	mod_atmos_phy_ch_vars.o      \
	mod_atmos_phy_rd_vars.o      \
	mod_atmos_phy_sf_vars.o      \
	mod_atmos_phy_tb_vars.o      \
	mod_atmos_phy_cp_vars.o      \
	mod_atmos_vars.o             \
	mod_atmos_dyn_driver.o       \
	mod_atmos_phy_mp_driver.o    \
	mod_atmos_phy_ae_driver.o    \
	mod_atmos_phy_ch_driver.o    \
	mod_atmos_phy_rd_driver.o    \
	mod_atmos_phy_sf_driver.o    \
	mod_atmos_phy_tb_driver.o    \
	mod_atmos_phy_cp_driver.o    \
	mod_atmos_driver.o           \
	mod_ocean_vars.o             \
	mod_ocean_phy_driver.o       \
	mod_ocean_driver.o           \
	mod_land_vars.o              \
	mod_land_phy_driver.o        \
	mod_land_driver.o            \
	mod_urban_vars.o             \
	mod_urban_phy_driver.o       \
	mod_urban_driver.o

OBJS_RUN =	\
	mod_les_driver.o

OBJS_INIT =	\
	mod_cnvtopo.o                \
	mod_copytopo.o               \
	mod_cnvlanduse.o             \
	mod_convert.o                \
	mod_mktopo.o                 \
	mod_realinput_scale.o        \
	mod_realinput_wrfarw.o       \
	mod_realinput_nicam.o        \
	mod_realinput_grads.o        \
	mod_realinput_scale.o        \
	mod_realinput.o              \
	mod_mkinit.o

all:
	$(MAKE) makedir
	$(MAKE) -C $(DCUTILSDIR)
	$(MAKE) -C $(GTOOLDIR)
	$(MAKE) -C $(SCALELIBDIR)/src
	@echo;echo "Entering scale-les...";echo "Current version is " $(VERSION)
	$(MAKE) makebin
	@echo "Complete making scale-les."

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

$(BINDIR)/$(BINNAME)  : $(BUILD_DIR)/$(BINNAME).o  $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(patsubst %,$(BUILD_DIR)/%,$(OBJS_RUN))  $(BUILD_DIR)/mod_user.o
	$(FC) $(LDFLAGS) -o $@ $^ $(DEPENDLIB) $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BINDIR)/$(INITNAME) : $(BUILD_DIR)/$(INITNAME).o $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(patsubst %,$(BUILD_DIR)/%,$(OBJS_INIT)) $(BUILD_DIR)/mod_user.o
	$(FC) $(LDFLAGS) -o $@ $^ $(DEPENDLIB) $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BINDIR)/$(PPNAME)   : $(BUILD_DIR)/$(PPNAME).o   $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(patsubst %,$(BUILD_DIR)/%,$(OBJS_INIT)) $(BUILD_DIR)/mod_user.o
	$(FC) $(LDFLAGS) -o $@ $^ $(DEPENDLIB) $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BUILD_DIR)/$(BINNAME).o  : $(BINNAME).f90  $(DEPENDLIB) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(patsubst %,$(BUILD_DIR)/%,$(OBJS_RUN))  $(BUILD_DIR)/mod_user.o
$(BUILD_DIR)/$(INITNAME).o : $(INITNAME).f90 $(DEPENDLIB) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(patsubst %,$(BUILD_DIR)/%,$(OBJS_INIT)) $(BUILD_DIR)/mod_user.o
$(BUILD_DIR)/$(PPNAME).o   : $(PPNAME).f90   $(DEPENDLIB) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(patsubst %,$(BUILD_DIR)/%,$(OBJS_INIT)) $(BUILD_DIR)/mod_user.o

$(BUILD_DIR)/mod_admin_restart.o          : mod_admin_restart.f90 $(DEPENDLIB) \
                                            mod_atmos_vars.o        \
                                            mod_atmos_dyn_vars.o    \
                                            mod_atmos_phy_ae_vars.o \
                                            mod_atmos_phy_ch_vars.o \
                                            mod_atmos_phy_cp_vars.o \
                                            mod_atmos_phy_mp_vars.o \
                                            mod_atmos_phy_rd_vars.o \
                                            mod_atmos_phy_sf_vars.o \
                                            mod_atmos_phy_tb_vars.o \
                                            mod_ocean_vars.o        \
                                            mod_land_vars.o         \
                                            mod_urban_vars.o
$(BUILD_DIR)/mod_admin_time.o             : mod_admin_time.f90 $(DEPENDLIB)

$(BUILD_DIR)/mod_atmos_admin.o            : mod_atmos_admin.f90 $(DEPENDLIB)

$(BUILD_DIR)/mod_atmos_dyn_vars.o         : mod_atmos_dyn_vars.f90    $(DEPENDLIB) mod_atmos_admin.o
$(BUILD_DIR)/mod_atmos_phy_cp_vars.o      : mod_atmos_phy_cp_vars.f90 $(DEPENDLIB) mod_atmos_admin.o
$(BUILD_DIR)/mod_atmos_phy_mp_vars.o      : mod_atmos_phy_mp_vars.f90 $(DEPENDLIB) mod_atmos_admin.o
$(BUILD_DIR)/mod_atmos_phy_rd_vars.o      : mod_atmos_phy_rd_vars.f90 $(DEPENDLIB) mod_atmos_admin.o
$(BUILD_DIR)/mod_atmos_phy_sf_vars.o      : mod_atmos_phy_sf_vars.f90 $(DEPENDLIB) mod_atmos_admin.o
$(BUILD_DIR)/mod_atmos_phy_tb_vars.o      : mod_atmos_phy_tb_vars.f90 $(DEPENDLIB) mod_atmos_admin.o
$(BUILD_DIR)/mod_atmos_phy_ch_vars.o      : mod_atmos_phy_ch_vars.f90 $(DEPENDLIB) mod_atmos_admin.o
$(BUILD_DIR)/mod_atmos_phy_ae_vars.o      : mod_atmos_phy_ae_vars.f90 $(DEPENDLIB) mod_atmos_admin.o
$(BUILD_DIR)/mod_atmos_vars.o             : mod_atmos_vars.f90        $(DEPENDLIB) mod_atmos_admin.o mod_atmos_dyn_vars.o mod_atmos_phy_mp_vars.o mod_atmos_phy_ae_vars.o mod_atmos_phy_ch_vars.o mod_atmos_phy_rd_vars.o mod_atmos_phy_sf_vars.o mod_atmos_phy_tb_vars.o mod_atmos_phy_cp_vars.o
$(BUILD_DIR)/mod_atmos_dyn_driver.o       : mod_atmos_dyn_driver.f90    $(DEPENDLIB) mod_atmos_admin.o mod_atmos_vars.o mod_atmos_dyn_vars.o
$(BUILD_DIR)/mod_atmos_phy_cp_driver.o    : mod_atmos_phy_cp_driver.f90 $(DEPENDLIB) mod_atmos_admin.o mod_atmos_vars.o mod_atmos_phy_cp_vars.o
$(BUILD_DIR)/mod_atmos_phy_mp_driver.o    : mod_atmos_phy_mp_driver.f90 $(DEPENDLIB) mod_atmos_admin.o mod_atmos_vars.o mod_atmos_phy_mp_vars.o
$(BUILD_DIR)/mod_atmos_phy_rd_driver.o    : mod_atmos_phy_rd_driver.f90 $(DEPENDLIB) mod_atmos_admin.o mod_atmos_vars.o mod_atmos_phy_rd_vars.o mod_cpl_admin.o mod_cpl_vars.o
$(BUILD_DIR)/mod_atmos_phy_sf_driver.o    : mod_atmos_phy_sf_driver.f90 $(DEPENDLIB) mod_atmos_admin.o mod_atmos_vars.o mod_atmos_phy_sf_vars.o mod_cpl_admin.o mod_cpl_vars.o
$(BUILD_DIR)/mod_atmos_phy_tb_driver.o    : mod_atmos_phy_tb_driver.f90 $(DEPENDLIB) mod_atmos_admin.o mod_atmos_vars.o mod_atmos_phy_tb_vars.o
$(BUILD_DIR)/mod_atmos_phy_ch_driver.o    : mod_atmos_phy_ch_driver.f90 $(DEPENDLIB) mod_atmos_admin.o mod_atmos_vars.o mod_atmos_phy_ch_vars.o
$(BUILD_DIR)/mod_atmos_phy_ae_driver.o    : mod_atmos_phy_ae_driver.f90 $(DEPENDLIB) mod_atmos_admin.o mod_atmos_vars.o mod_atmos_phy_ae_vars.o
$(BUILD_DIR)/mod_atmos_driver.o           : mod_atmos_driver.f90        $(DEPENDLIB) mod_atmos_admin.o mod_atmos_vars.o mod_admin_time.o \
                                            mod_atmos_dyn_driver.o    \
                                            mod_atmos_phy_cp_driver.o \
                                            mod_atmos_phy_mp_driver.o \
                                            mod_atmos_phy_rd_driver.o \
                                            mod_atmos_phy_sf_driver.o \
                                            mod_atmos_phy_tb_driver.o \
                                            mod_atmos_phy_ch_driver.o \
                                            mod_atmos_phy_ae_driver.o

$(BUILD_DIR)/mod_ocean_admin.o            : mod_ocean_admin.f90      $(DEPENDLIB)
$(BUILD_DIR)/mod_ocean_vars.o             : mod_ocean_vars.f90       $(DEPENDLIB) mod_ocean_admin.o
$(BUILD_DIR)/mod_ocean_phy_driver.o       : mod_ocean_phy_driver.f90 $(DEPENDLIB) mod_ocean_admin.o mod_ocean_vars.o mod_admin_restart.o
$(BUILD_DIR)/mod_ocean_driver.o           : mod_ocean_driver.f90     $(DEPENDLIB) mod_ocean_admin.o mod_ocean_vars.o mod_cpl_vars.o \
                                            mod_ocean_phy_driver.o

$(BUILD_DIR)/mod_land_admin.o             : mod_land_admin.f90      $(DEPENDLIB)
$(BUILD_DIR)/mod_land_vars.o              : mod_land_vars.f90       $(DEPENDLIB) mod_land_admin.o
$(BUILD_DIR)/mod_land_phy_driver.o        : mod_land_phy_driver.f90 $(DEPENDLIB) mod_land_admin.o mod_land_vars.o mod_admin_restart.o
$(BUILD_DIR)/mod_land_driver.o            : mod_land_driver.f90     $(DEPENDLIB) mod_land_admin.o mod_land_vars.o mod_cpl_vars.o \
                                            mod_land_phy_driver.o

$(BUILD_DIR)/mod_urban_admin.o            : mod_urban_admin.f90      $(DEPENDLIB)
$(BUILD_DIR)/mod_urban_vars.o             : mod_urban_vars.f90       $(DEPENDLIB) mod_urban_admin.o
$(BUILD_DIR)/mod_urban_phy_driver.o       : mod_urban_phy_driver.f90 $(DEPENDLIB) mod_urban_admin.o mod_urban_vars.o mod_admin_restart.o
$(BUILD_DIR)/mod_urban_driver.o           : mod_urban_driver.f90     $(DEPENDLIB) mod_urban_admin.o mod_urban_vars.o mod_cpl_vars.o \
                                            mod_urban_phy_driver.o

$(BUILD_DIR)/mod_cpl_admin.o              : mod_cpl_admin.f90  $(DEPENDLIB) mod_atmos_admin.o mod_ocean_admin.o mod_land_admin.o mod_urban_admin.o
$(BUILD_DIR)/mod_cpl_vars.o               : mod_cpl_vars.f90   $(DEPENDLIB) mod_cpl_admin.f90
$(BUILD_DIR)/mod_cpl_driver.o             : mod_cpl_driver.f90 $(DEPENDLIB) mod_cpl_admin.o mod_cpl_vars.o

$(BUILD_DIR)/mod_mktopo.o                 : mod_mktopo.f90     $(DEPENDLIB)
$(BUILD_DIR)/mod_mkinit.o                 : mod_mkinit.f90     $(DEPENDLIB) mod_atmos_vars.o mod_ocean_vars.o mod_land_vars.o mod_urban_vars.o mod_realinput.o
$(BUILD_DIR)/mod_realinput_scale.o        : mod_realinput_scale.f90 $(DEPENDLIB)
$(BUILD_DIR)/mod_realinput_wrfarw.o       : mod_realinput_wrfarw.f90 $(DEPENDLIB)
$(BUILD_DIR)/mod_realinput_nicam.o        : mod_realinput_nicam.f90 $(DEPENDLIB)
$(BUILD_DIR)/mod_realinput_grads.o        : mod_realinput_grads.f90 $(DEPENDLIB)
$(BUILD_DIR)/mod_realinput.o              : mod_realinput.f90  $(DEPENDLIB) mod_atmos_vars.o mod_ocean_vars.o mod_land_vars.o mod_urban_vars.o mod_realinput_scale.o mod_realinput_wrfarw.o mod_realinput_nicam.o mod_realinput_grads.o
$(BUILD_DIR)/mod_convert.o                : mod_convert.f90    $(DEPENDLIB) mod_cnvtopo.o mod_cnvlanduse.o
$(BUILD_DIR)/mod_cnvtopo.o                : mod_cnvtopo.f90    $(DEPENDLIB) mod_copytopo.o
$(BUILD_DIR)/mod_copytopo.o               : mod_copytopo.f90   $(DEPENDLIB)
$(BUILD_DIR)/mod_cnvlanduse.o             : mod_cnvlanduse.f90 $(DEPENDLIB)

$(BUILD_DIR)/mod_user.o                   : mod_user.f90 $(DEPENDLIB) $(patsubst %,$(BUILD_DIR)/%,$(OBJS))

$(BUILD_DIR)/mod_les_driver.o             : mod_les_driver.f90 $(DEPENDLIB) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o


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



.SUFFIXES:
.SUFFIXES: .o .f90 .c .mod

$(BUILD_DIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MODDIR) -I$(SCALELIBDIR)/include -I$(SCALELESDIR)/include $(PAPI_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)

%.mod: %.f90
	make $(patsubst %.f90,%.o,$<)


.PHONY: clean dcutils gtool scalelib
