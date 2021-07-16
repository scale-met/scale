################################################################################
#
# Makefile for main program
#
################################################################################

TOPDIR     = $(abspath ../..)
SYSDEP_DIR = $(abspath ../../sysdep)

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude


BINNAME  = scale-rm
INITNAME = scale-rm_init
PPNAME   = scale-rm_pp

BINS = $(BINDIR)/$(BINNAME)$(POSTFIX)  \
       $(BINDIR)/$(INITNAME)$(POSTFIX) \
       $(BINDIR)/$(PPNAME)$(POSTFIX)
LIBS = $(LIBDIR)/libscale$(POSTFIX).a

ifneq ($(SCALE_IGNORE_SRCDEP),T)
  DEPENDLIB = $(LIBS)
endif

VPATH = $(BUILD_DIR): \
        admin:        \
        atmos:        \
        ocean:        \
        land:         \
        urban:        \
        lake:         \
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
	mod_admin_versioncheck.o  \
	mod_admin_time.o          \
	mod_atmos_admin.o         \
	mod_ocean_admin.o         \
	mod_land_admin.o          \
	mod_urban_admin.o         \
	mod_lake_admin.o          \
	mod_cpl_admin.o           \
	\
	mod_atmos_dyn_vars.o      \
	mod_atmos_phy_cp_vars.o   \
	mod_atmos_phy_lt_vars.o   \
	mod_atmos_phy_mp_vars.o   \
	mod_atmos_phy_rd_vars.o   \
	mod_atmos_phy_sf_vars.o   \
	mod_atmos_phy_tb_vars.o   \
	mod_atmos_phy_bl_vars.o   \
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
	mod_atmos_phy_lt_driver.o \
	mod_atmos_phy_mp_driver.o \
	mod_atmos_phy_ae_driver.o \
	mod_atmos_phy_ch_driver.o \
	mod_atmos_phy_rd_driver.o \
	mod_atmos_phy_sf_driver.o \
	mod_atmos_phy_tb_driver.o \
	mod_atmos_phy_bl_driver.o \
	mod_atmos_phy_cp_driver.o \
	mod_atmos_bnd_driver.o	  \
	mod_atmos_driver.o        \
	mod_ocean_driver.o        \
	mod_land_driver.o         \
	mod_urban_driver.o        \
	mod_cpl_driver.o          \
	\
	mod_copytopo.o            \
	mod_cnv2d.o               \
	mod_cnvtopo.o             \
	mod_cnvlanduse.o          \
	mod_cnvuser.o             \
	mod_convert.o             \
	mod_mktopo.o              \
	mod_realinput_scale.o     \
	mod_realinput_wrfarw.o    \
	mod_realinput_grads.o     \
	mod_realinput.o           \
	mod_mkinit.o
#	mod_realinput_nicam.o     \

MAIN_OBJS =	\
	mod_launcher.o            \
	mod_rm_driver.o           \
	mod_rm_prep.o

OBJS = 	\
	$(SUB_OBJS)	\
	$(MAIN_OBJS)

all:
	$(MAKE) build
	$(MAKE) install

build:
	$(MAKE) info
	$(MAKE) makelib
	@echo;echo "Entering scale-rm...";echo "Current version is " $(VERSION)
	mkdir -p $(BUILD_DIR)
	$(MAKE) makebin
	@echo "Complete making scale-rm."

install:
	mkdir -p $(BINDIR)
	$(MAKE) $(BINS)

info:
	@$(MAKE) -C $(SCALELIBDIR)/src --no-print-directory info
	@echo
	@echo "SCALE-RM"
	@$(MAKE) --no-print-directory conflog

makelib:
	$(MAKE) -C $(SCALELIBDIR)/src

makebin: $(BUILD_DIR)/$(BINNAME)$(POSTFIX) $(BUILD_DIR)/$(PPNAME)$(POSTFIX)

$(BUILD_DIR)/$(BINNAME)$(POSTFIX)  : $(BUILD_DIR)/$(BINNAME).o $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $< $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o -L$(LIBDIR) -lscale$(POSTFIX) $(CONTRIB_LIBS) $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BUILD_DIR)/$(PPNAME)$(POSTFIX) : $(BUILD_DIR)/$(PPNAME).o $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $< $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o -L$(LIBDIR) -lscale$(POSTFIX) $(CONTRIB_LIBS) $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)

$(BINDIR)/$(BINNAME)$(POSTFIX)  : $(BUILD_DIR)/$(BINNAME)$(POSTFIX)
	install $< $@

$(BINDIR)/$(PPNAME)$(POSTFIX)   : $(BUILD_DIR)/$(PPNAME)$(POSTFIX)
	install $< $@

$(BINDIR)/$(INITNAME)$(POSTFIX) : $(BINDIR)/$(PPNAME)$(POSTFIX)
	(cd $(BINDIR); ln -sf `basename $<` $@)



ifneq ($(SCALE_IGNORE_SRCDEP),T)
  $(BUILD_DIR)/$(BINNAME).o : $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o
  $(BUILD_DIR)/$(PPNAME).o  : $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o
  $(BUILD_DIR)/mod_user.o   : $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(SUB_OBJS))
endif


allclean: distclean
	$(MAKE) -C $(SCALELIBDIR)/src allclean
	rm -f $(BINDIR)/$(BINNAME)*
	rm -f $(BINDIR)/$(INITNAME)*
	rm -f $(BINDIR)/$(PPNAME)*
	rm -rf $(TOPDIR)/bin


distclean: clean
	rm -f $(BINDIR)/$(BINNAME)$(POSTFIX)
	rm -f $(BINDIR)/$(INITNAME)$(POSTFIX)
	rm -f $(BINDIR)/$(PPNAME)$(POSTFIX)
	rm -rf ./.libs*

clean:
	rm -rf $(BUILD_DIR)
	rm -f *.o *.mod *.lst *~

dependency:
	./makedepend .



.SUFFIXES:
.SUFFIXES: .o .F90 .c .mod

%.mod: %.F90
	make $(patsubst %.F90,%.o,$<)

$(BUILD_DIR)/%.o : %.F90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(MODDIR) $(CONTRIB_INCLUDE) -I$(SCALELIBDIR)/include -I$(SCALERMDIR)/include $(SCALE_PAPI_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

$(BUILD_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONY : clean distclean allclean scalelib depend dependency

include $(TOPDIR)/utils/make/Make_environments

ifneq ($(SCALE_IGNORE_SRCDEP),T)
-include depend
endif
