################################################################################
#
# Makefile for main program
#
################################################################################

TOPDIR     = $(abspath ../..)
BUILD_DIR  = ./.libs
SYSDEP_DIR = ../../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

BINNAME  = scale-gm
PPNAME1  = gm_mkmnginfo
PPNAME2  = gm_mkrawgrid
PPNAME3  = gm_mkhgrid
PPNAME4  = gm_mkvlayer
PPNAME5  = gm_mkllmap
PPNAME6  = gm_fio_dump
PPNAME7  = gm_fio_cat
PPNAME8  = gm_fio_sel
PPNAME9  = gm_fio_ico2ll

BINS = $(BINDIR)/$(BINNAME) \
       $(BINDIR)/$(PPNAME1) \
       $(BINDIR)/$(PPNAME2) \
       $(BINDIR)/$(PPNAME3) \
       $(BINDIR)/$(PPNAME4) \
       $(BINDIR)/$(PPNAME5) \
       $(BINDIR)/$(PPNAME6) \
       $(BINDIR)/$(PPNAME7) \
       $(BINDIR)/$(PPNAME8) \
       $(BINDIR)/$(PPNAME9)
LIBS = $(LIBDIR)/libscale.a

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
	mod_adm.o            \
	mod_fio.o            \
	mod_comm.o           \
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
	mod_mnginfo_light.o  \
	mod_option.o         \
	mod_netcdf.o

all:
	$(MAKE) envlog
	$(MAKE) makedir
	$(MAKE) -C $(SCALELIBDIR)/src
	@echo;echo "Entering scale-gm...";echo "Current version is " $(VERSION)
	$(MAKE) makebin
	@echo "Complete making scale-gm."

makedir:
	mkdir -p $(BINDIR)
	mkdir -p $(LIBDIR)
	mkdir -p $(BUILD_DIR)

$(LIBDIR)/libscale.a:
	$(MAKE) -C $(SCALELIBDIR)/src

makebin: $(BINS)

$(BINDIR)/$(BINNAME) : $(BUILD_DIR)/scale-gm.o       $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BINDIR)/$(PPNAME1) : $(BUILD_DIR)/prg_mkmnginfo.o  $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BINDIR)/$(PPNAME2) : $(BUILD_DIR)/prg_mkrawgrid.o  $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BINDIR)/$(PPNAME3) : $(BUILD_DIR)/prg_mkhgrid.o    $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BINDIR)/$(PPNAME4) : $(BUILD_DIR)/prg_mkvlayer.o   $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BINDIR)/$(PPNAME5) : $(BUILD_DIR)/prg_mkllmap.o    $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BINDIR)/$(PPNAME6) : $(BUILD_DIR)/prg_fio_dump.o   $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BINDIR)/$(PPNAME7) : $(BUILD_DIR)/prg_fio_cat.o    $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BINDIR)/$(PPNAME8) : $(BUILD_DIR)/prg_fio_sel.o    $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)
$(BINDIR)/$(PPNAME9) : $(BUILD_DIR)/prg_fio_ico2ll.o $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)

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



allclean: distclean
	$(MAKE) -C $(SCALELIBDIR)/src allclean
	rm -f $(BINDIR)/$(BINNAME)
	rm -f $(BINDIR)/$(PPNAME1)
	rm -f $(BINDIR)/$(PPNAME2)
	rm -f $(BINDIR)/$(PPNAME3)
	rm -f $(BINDIR)/$(PPNAME4)
	rm -f $(BINDIR)/$(PPNAME5)
	rm -f $(BINDIR)/$(PPNAME6)
	rm -f $(BINDIR)/$(PPNAME7)
	rm -f $(BINDIR)/$(PPNAME8)
	rm -f $(BINDIR)/$(PPNAME9)

distclean: clean
	rm -f $(BUILD_DIR)/$(BINNAME)
	rm -f $(BUILD_DIR)/$(PPNAME1)
	rm -f $(BUILD_DIR)/$(PPNAME2)
	rm -f $(BUILD_DIR)/$(PPNAME3)
	rm -f $(BUILD_DIR)/$(PPNAME4)
	rm -f $(BUILD_DIR)/$(PPNAME5)
	rm -f $(BUILD_DIR)/$(PPNAME6)
	rm -f $(BUILD_DIR)/$(PPNAME7)
	rm -f $(BUILD_DIR)/$(PPNAME8)
	rm -f $(BUILD_DIR)/$(PPNAME9)

clean:
	rm -f $(BUILD_DIR)/*.*

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
