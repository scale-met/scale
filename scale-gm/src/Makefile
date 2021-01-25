################################################################################
#
# Makefile for main program
#
################################################################################

TOPDIR     = $(abspath ../..)
SYSDEP_DIR = $(abspath ../../sysdep)

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

BINNAME  = scale-gm
PPNAME01 = gm_mkmnginfo
PPNAME02 = gm_mkrawgrid
PPNAME03 = gm_mkhgrid
PPNAME04 = gm_mkvlayer
PPNAME05 = gm_mkllmap
PPNAME06 = gm_fio_dump
PPNAME07 = gm_fio_cat
PPNAME08 = gm_fio_sel
PPNAME09 = gm_fio_ico2ll
PPNAME10 = gm_ll2ico

BINS = $(BINDIR)/$(BINNAME)$(POSTFIX)  \
       $(BINDIR)/$(PPNAME01)$(POSTFIX) \
       $(BINDIR)/$(PPNAME02)$(POSTFIX) \
       $(BINDIR)/$(PPNAME03)$(POSTFIX) \
       $(BINDIR)/$(PPNAME04)$(POSTFIX) \
       $(BINDIR)/$(PPNAME05)$(POSTFIX) \
       $(BINDIR)/$(PPNAME06)$(POSTFIX) \
       $(BINDIR)/$(PPNAME07)$(POSTFIX) \
       $(BINDIR)/$(PPNAME08)$(POSTFIX) \
       $(BINDIR)/$(PPNAME09)$(POSTFIX) \
       $(BINDIR)/$(PPNAME10)$(POSTFIX)
LIBS = $(LIBDIR)/libscale$(POSTFIX).a

ifneq ($(SCALE_IGNORE_SRCDEP),T)
  DEPENDLIB = $(LIBS)
endif

VPATH = $(BUILD_DIR):           \
        admin:                  \
        atmos:                  \
        atmos/dcmip:            \
        atmos/forcing:          \
        atmos/share:            \
        atmos-gm/communication: \
        atmos-gm/dynamics:      \
        atmos-gm/grid:          \
        common:                 \
        coupler:                \
        io:                     \
        ocean:                  \
        land:                   \
        lake:                   \
        urban:                  \
        preprocess:             \
        tool:


VERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(VERSION),)
  VERSION  = $(shell cat VERSION)
else
  VERSION := $(VERSION)
endif

OBJS =	\
	fio.o                \
	fiof.o               \
	mod_io_param.o       \
	mod_fio.o            \
	mod_ideal_topo.o     \
	mod_grd.o            \
	mod_gmtr.o           \
	mod_oprt.o           \
	mod_vmtr.o           \
	mod_time.o           \
	dcmip_initial_conditions_test_1_2_3_v5.o \
	baroclinic_wave_test.o                   \
	supercell_test.o                         \
	tropical_cyclone_test.o                  \
	simple_physics_v6_1.o                    \
	kessler.o                                \
	Terminator.o                             \
	mod_chemvar.o        \
	mod_ocean_admin.o   \
	mod_land_admin.o   \
	mod_lake_admin.o   \
	mod_urban_admin.o   \
	mod_atmos_admin.o    \
	mod_cpl_vars.o      \
	mod_cpl_admin.o     \
	mod_runconf.o        \
	mod_gtl.o            \
	mod_gm_statistics.o  \
	mod_vintrpl.o        \
	mod_bndcnd.o         \
	mod_cnvvar.o         \
	mod_src.o            \
	mod_ideal_init.o     \
	mod_prgvar.o         \
	mod_bsstate.o        \
	mod_history.o        \
	mod_history_vars.o   \
	mod_embudget.o       \
	mod_extdata.o        \
	mod_af_heldsuarez.o  \
	mod_grd_conversion.o \
	mod_atmos_vars.o     \
	mod_atmos_phy_mp_vars.o     \
	mod_atmos_phy_ae_vars.o     \
	mod_atmos_phy_rd_vars.o     \
	mod_atmos_phy_sf_vars.o     \
	mod_atmos_phy_bl_vars.o     \
	mod_mp_vars.o	     \
	mod_af_trcadv.o      \
	mod_oprt3d.o         \
	mod_numfilter.o      \
	mod_vi.o             \
	mod_src_tracer.o     \
	mod_nudge.o          \
	mod_gm_topography.o  \
	mod_atmos_phy_driver.o	\
	mod_atmos_phy_mp_driver.o   \
	mod_atmos_phy_rd_driver.o \
	mod_atmos_phy_sf_driver.o \
	mod_atmos_phy_bl_driver.o \
	mod_atmos_surface.o     \
	mod_dynamics.o       \
	mod_af_dcmip.o       \
	mod_forcing_driver.o \
	mod_gm_driver.o      \
	orgico.o             \
	mod_latlon.o         \
	mod_mkgrd.o          \
	mod_gm_convert.o     \
	mod_gm_cnvtopo.o     \
	mod_gm_cnv2d.o       \
	mod_mnginfo_light.o  \
	mod_option.o         \
	mod_netcdf.o

all:
	$(MAKE) build
	$(MAKE) install

build:
	$(MAKE) info
	$(MAKE) makelib
	@echo;echo "Entering scale-gm...";echo "Current version is " $(VERSION)
	mkdir -p $(BUILD_DIR)
	$(MAKE) $(subst $(BINDIR),$(BUILD_DIR),$(BINS))
	@echo "Complete making scale-gm."

install:
	mkdir -p $(BINDIR)
	$(MAKE) $(BINS)

makelib:
	$(MAKE) -C $(SCALELIBDIR)/src

info:
	@$(MAKE) -C $(SCALELIBDIR)/src --no-print-directory info
	@echo
	@echo "SCALE-GM"
	@$(MAKE) --no-print-directory conflog

$(BUILD_DIR)/$(BINNAME)$(POSTFIX) : $(BUILD_DIR)/scale-gm.o        $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BUILD_DIR)/$(PPNAME01)$(POSTFIX) : $(BUILD_DIR)/prg_mkmnginfo.o  $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BUILD_DIR)/$(PPNAME02)$(POSTFIX) : $(BUILD_DIR)/prg_mkrawgrid.o  $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BUILD_DIR)/$(PPNAME03)$(POSTFIX) : $(BUILD_DIR)/prg_mkhgrid.o    $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BUILD_DIR)/$(PPNAME04)$(POSTFIX) : $(BUILD_DIR)/prg_mkvlayer.o   $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BUILD_DIR)/$(PPNAME05)$(POSTFIX) : $(BUILD_DIR)/prg_mkllmap.o    $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BUILD_DIR)/$(PPNAME06)$(POSTFIX) : $(BUILD_DIR)/prg_fio_dump.o   $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BUILD_DIR)/$(PPNAME07)$(POSTFIX) : $(BUILD_DIR)/prg_fio_cat.o    $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BUILD_DIR)/$(PPNAME08)$(POSTFIX) : $(BUILD_DIR)/prg_fio_sel.o    $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BUILD_DIR)/$(PPNAME09)$(POSTFIX) : $(BUILD_DIR)/prg_fio_ico2ll.o $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BUILD_DIR)/$(PPNAME10)$(POSTFIX) : $(BUILD_DIR)/prg_ll2ico.o     $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)

$(BINDIR)/$(BINNAME)$(POSTFIX) : $(BUILD_DIR)/$(BINNAME)$(POSTFIX)
	install $< $@
$(BINDIR)/$(PPNAME01)$(POSTFIX) : $(BUILD_DIR)/$(PPNAME01)$(POSTFIX)
	install $< $@
$(BINDIR)/$(PPNAME02)$(POSTFIX) : $(BUILD_DIR)/$(PPNAME02)$(POSTFIX)
	install $< $@
$(BINDIR)/$(PPNAME03)$(POSTFIX) : $(BUILD_DIR)/$(PPNAME03)$(POSTFIX)
	install $< $@
$(BINDIR)/$(PPNAME04)$(POSTFIX) : $(BUILD_DIR)/$(PPNAME04)$(POSTFIX)
	install $< $@
$(BINDIR)/$(PPNAME05)$(POSTFIX) : $(BUILD_DIR)/$(PPNAME05)$(POSTFIX)
	install $< $@
$(BINDIR)/$(PPNAME06)$(POSTFIX) : $(BUILD_DIR)/$(PPNAME06)$(POSTFIX)
	install $< $@
$(BINDIR)/$(PPNAME07)$(POSTFIX) : $(BUILD_DIR)/$(PPNAME07)$(POSTFIX)
	install $< $@
$(BINDIR)/$(PPNAME08)$(POSTFIX) : $(BUILD_DIR)/$(PPNAME08)$(POSTFIX)
	install $< $@
$(BINDIR)/$(PPNAME09)$(POSTFIX) : $(BUILD_DIR)/$(PPNAME09)$(POSTFIX)
	install $< $@
$(BINDIR)/$(PPNAME10)$(POSTFIX) : $(BUILD_DIR)/$(PPNAME10)$(POSTFIX)
	install $< $@


$(BUILD_DIR)/scale-gm.o       : scale-gm.F90       $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_mkmnginfo.o  : prg_mkmnginfo.F90  $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_mkrawgrid.o  : prg_mkrawgrid.F90  $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_mkhgrid.o    : prg_mkhgrid.F90    $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_mkvlayer.o   : prg_mkvlayer.F90   $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_mkllmap.o    : prg_mkllmap.F90    $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_fio_dump.o   : prg_fio_dump.F90   $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_fio_cat.o    : prg_fio_cat.F90    $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_fio_sel.o    : prg_fio_sel.F90    $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_fio_ico2ll.o : prg_fio_ico2ll.F90 $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_ll2ico.o     : prg_ll2ico.F90     $(patsubst %,$(BUILD_DIR)/%,$(OBJS))



allclean: distclean
	$(MAKE) -C $(SCALELIBDIR)/src allclean
	rm -f $(BINDIR)/$(BINNAME)*
	rm -f $(BINDIR)/$(PPNAME01)*
	rm -f $(BINDIR)/$(PPNAME02)*
	rm -f $(BINDIR)/$(PPNAME03)*
	rm -f $(BINDIR)/$(PPNAME04)*
	rm -f $(BINDIR)/$(PPNAME05)*
	rm -f $(BINDIR)/$(PPNAME06)*
	rm -f $(BINDIR)/$(PPNAME07)*
	rm -f $(BINDIR)/$(PPNAME08)*
	rm -f $(BINDIR)/$(PPNAME09)*
	rm -f $(BINDIR)/$(PPNAME10)*


distclean: clean
	rm -f $(BINDIR)/$(BINNAME)$(POSTFIX)
	rm -f $(BINDIR)/$(PPNAME01)$(POSTFIX)
	rm -f $(BINDIR)/$(PPNAME02)$(POSTFIX)
	rm -f $(BINDIR)/$(PPNAME03)$(POSTFIX)
	rm -f $(BINDIR)/$(PPNAME04)$(POSTFIX)
	rm -f $(BINDIR)/$(PPNAME05)$(POSTFIX)
	rm -f $(BINDIR)/$(PPNAME06)$(POSTFIX)
	rm -f $(BINDIR)/$(PPNAME07)$(POSTFIX)
	rm -f $(BINDIR)/$(PPNAME08)$(POSTFIX)
	rm -f $(BINDIR)/$(PPNAME09)$(POSTFIX)
	rm -f $(BINDIR)/$(PPNAME10)$(POSTFIX)
	rm -rf ./libs*

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
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(MODDIR) $(CONTRIB_INCLUDE) -I$(SCALELIBDIR)/include -I$(SCALEGMDIR)/include $(SCALE_PAPI_INCLUDE) $(SCALE_NETCDF_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<
$(BUILD_DIR)/mod_oprt.o :
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) $(FFLAGS_TENTATIVE) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(MODDIR) $(CONTRIB_INCLUDE) -I$(SCALELIBDIR)/include -I$(SCALEGMDIR)/include $(SCALE_PAPI_INCLUDE) $(SCALE_NETCDF_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

$(BUILD_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -I$(SCALEGMDIR)/include -o $@ -c $<

.PHONY : clean distclean allclean scalelib depend dependency

include $(TOPDIR)/utils/make/Make_environments

ifneq ($(SCALE_IGNORE_SRCDEP),T)
-include depend
endif
