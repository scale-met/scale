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

BINNAME = scale-gm
PPNAME1 = gm_mkmnginfo
PPNAME2 = gm_mkrawgrid
PPNAME3 = gm_mkhgrid
PPNAME4 = gm_mkvlayer
PPNAME5 = gm_mkllmap
PPNAME6 = gm_fio_dump
PPNAME7 = gm_fio_cat
PPNAME8 = gm_fio_sel
PPNAME9 = gm_fio_ico2ll

#ifeq ($(NOINIT),)
#  BINS = $(BINDIR)/$(BINNAME) $(BINDIR)/$(PPNAME1) $(BINDIR)/$(PPNAME2) $(BINDIR)/$(PPNAME3) $(BINDIR)/$(PPNAME4) $(BINDIR)/$(PPNAME5) $(BINDIR)/$(PPNAME6) $(BINDIR)/$(PPNAME7) $(BINDIR)/$(PPNAME8) $(BINDIR)/$(PPNAME9)
#else
  BINS = $(BINDIR)/$(BINNAME)
#endif

DEPENDLIB = $(LIBDIR)/libscale.a $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a

VPATH = $(BUILD_DIR):           \
        admin:                  \
        atmos:                  \
        atmos/dcmip:            \
        atmos-gm/communication: \
        atmos-gm/dynamics:      \
        atmos-gm/grid:          \
        atmos-physics/forcing:  \
        atmos-physics/share:    \
        common:                 \
        io:                     \
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
	mod_chemvar.o        \
	mod_runconf.o        \
	mod_gtl.o            \
	mod_gm_statistics.o  \
	mod_vintrpl.o        \
	mod_thrmdyn.o        \
	mod_bndcnd.o         \
	mod_cnvvar.o         \
	mod_src.o            \
	mod_ideal_init.o     \
	mod_prgvar.o         \
	mod_history.o        \
	mod_bsstate.o        \
	mod_history_vars.o   \
	mod_embudget.o       \
	mod_extdata.o        \
	mod_af_heldsuarez.o  \
	mod_af_trcadv.o      \
	mod_forcing_driver.o \
	mod_oprt3d.o         \
	mod_numfilter.o      \
	mod_vi.o             \
	mod_src_tracer.o     \
	mod_nudge.o          \
	mod_dynamics.o       \
	mod_gm_driver.o      \
	orgico.o             \
	mod_latlon.o         \
	mod_mkgrd.o          \
	mod_mnginfo_light.o  \
	mod_option.o         \
	mod_netcdf.o



all:
	$(MAKE) makedir
	$(MAKE) -C $(DCUTILSDIR)
	$(MAKE) -C $(GTOOLDIR)
	$(MAKE) -C $(SCALELIBDIR)/src
	@echo;echo "Entering scale-gm...";echo "Current version is " $(VERSION)
	$(MAKE) makebin
	@echo "Complete making scale-gm."

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

$(BINDIR)/$(BINNAME) : $(BUILD_DIR)/scale-gm.o       $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^ $(DEPENDLIB) $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BINDIR)/$(PPNAME1) : $(BUILD_DIR)/prg_mkmnginfo.o  $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^ $(DEPENDLIB) $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BINDIR)/$(PPNAME2) : $(BUILD_DIR)/prg_mkrawgrid.o  $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^ $(DEPENDLIB) $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BINDIR)/$(PPNAME3) : $(BUILD_DIR)/prg_mkhgrid.o    $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^ $(DEPENDLIB) $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BINDIR)/$(PPNAME4) : $(BUILD_DIR)/prg_mkvlayer.o   $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^ $(DEPENDLIB) $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BINDIR)/$(PPNAME5) : $(BUILD_DIR)/prg_mkllmap.o    $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^ $(DEPENDLIB) $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BINDIR)/$(PPNAME6) : $(BUILD_DIR)/prg_fio_dump.o   $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^ $(DEPENDLIB) $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BINDIR)/$(PPNAME7) : $(BUILD_DIR)/prg_fio_cat.o    $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^ $(DEPENDLIB) $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BINDIR)/$(PPNAME8) : $(BUILD_DIR)/prg_fio_sel.o    $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^ $(DEPENDLIB) $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BINDIR)/$(PPNAME9) : $(BUILD_DIR)/prg_fio_ico2ll.o $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^ $(DEPENDLIB) $(NETCDF_LIBS) $(LAPACK_LIBS) $(PAPI_LIBS)

$(BUILD_DIR)/scale-gm.o       : scale-gm.f90       $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_mkmnginfo.o  : prg_mkmnginfo.f90  $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_mkrawgrid.o  : prg_mkrawgrid.f90  $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_mkhgrid.o    : prg_mkhgrid.f90    $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_mkvlayer.o   : prg_mkvlayer.f90   $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_mkllmap.o    : prg_mkllmap.f90    $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_fio_dump.o   : prg_fio_dump.f90   $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_fio_cat.o    : prg_fio_cat.f90    $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_fio_sel.o    : prg_fio_sel.f90    $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
$(BUILD_DIR)/prg_fio_ico2ll.o : prg_fio_ico2ll.f90 $(patsubst %,$(BUILD_DIR)/%,$(OBJS))

$(BUILD_DIR)/fio.o                : fio.c $(SCALEGMDIR)/include/fio.h $(SCALEGMDIR)/include/fio_def.h
$(BUILD_DIR)/fiof.o               : fiof.c $(SCALEGMDIR)/include/fiof.h $(SCALEGMDIR)/include/fio.h $(SCALEGMDIR)/include/fio_def.h
$(BUILD_DIR)/orgico.o             : orgico.c $(SCALEGMDIR)/include/orgico.h $(SCALEGMDIR)/include/fio_def.h
$(BUILD_DIR)/mod_adm.o            : mod_adm.f90            mod_misc.o
$(BUILD_DIR)/mod_fio.o            : mod_fio.f90            mod_adm.o fio.o fiof.o
$(BUILD_DIR)/mod_comm.o           : mod_comm.f90           mod_adm.o
$(BUILD_DIR)/mod_grd.o            : mod_grd.f90            mod_adm.o mod_comm.o mod_fio.o mod_ideal_topo.o
$(BUILD_DIR)/mod_ideal_topo.o     : mod_ideal_topo.f90     mod_adm.o
$(BUILD_DIR)/mod_gmtr.o           : mod_gmtr.f90           mod_adm.o mod_comm.o mod_grd.o
$(BUILD_DIR)/mod_oprt.o           : mod_oprt.f90           mod_adm.o mod_grd.o mod_gmtr.o mod_comm.o
$(BUILD_DIR)/mod_oprt3d.o         : mod_oprt3d.f90         mod_adm.o mod_grd.o mod_gmtr.o mod_vmtr.o
$(BUILD_DIR)/mod_vmtr.o           : mod_vmtr.f90           mod_adm.o mod_grd.o mod_comm.o mod_gmtr.o mod_oprt.o
$(BUILD_DIR)/mod_time.o           : mod_time.f90           mod_adm.o
$(BUILD_DIR)/mod_gtl.o            : mod_gtl.f90            mod_adm.o
$(BUILD_DIR)/mod_gm_statistics.o  : mod_gm_statistics.f90  mod_adm.o mod_comm.o mod_gmtr.o mod_vmtr.o
$(BUILD_DIR)/mod_vintrpl.o        : mod_vintrpl.f90        mod_adm.o mod_grd.o
$(BUILD_DIR)/mod_latlon.o         : mod_latlon.f90         mod_adm.o mod_grd.o mod_oprt.o mod_comm.o mod_gmtr.o
$(BUILD_DIR)/mod_runconf.o        : mod_runconf.f90        mod_adm.o mod_chemvar.o
$(BUILD_DIR)/mod_cnvvar.o         : mod_cnvvar.f90         mod_adm.o mod_grd.o mod_thrmdyn.o mod_runconf.o mod_bndcnd.o mod_vmtr.o
$(BUILD_DIR)/mod_thrmdyn.o        : mod_thrmdyn.f90        mod_adm.o mod_runconf.o
$(BUILD_DIR)/mod_bndcnd.o         : mod_bndcnd.f90         mod_adm.o mod_grd.o mod_runconf.o
$(BUILD_DIR)/mod_bsstate.o        : mod_bsstate.f90        mod_adm.o mod_grd.o mod_runconf.o mod_vintrpl.o mod_bndcnd.o mod_thrmdyn.o
$(BUILD_DIR)/mod_history.o        : mod_history.f90        mod_adm.o mod_time.o mod_runconf.o mod_grd.o mod_fio.o mod_comm.o mod_gm_statistics.o mod_vintrpl.o mod_cnvvar.o mod_thrmdyn.o mod_vmtr.o mod_prgvar.o
$(BUILD_DIR)/mod_prgvar.o         : mod_prgvar.f90         mod_adm.o mod_runconf.o mod_comm.o mod_ideal_init.o mod_gm_statistics.o
$(BUILD_DIR)/mod_chemvar.o        : mod_chemvar.f90        mod_adm.o
$(BUILD_DIR)/mod_extdata.o        : mod_extdata.f90        mod_adm.o mod_fio.o mod_gtl.o mod_comm.o
$(BUILD_DIR)/mod_embudget.o       : mod_embudget.f90       mod_adm.o mod_time.o mod_grd.o mod_vmtr.o mod_gm_statistics.o mod_runconf.o mod_cnvvar.o mod_bsstate.o mod_thrmdyn.o mod_prgvar.o
$(BUILD_DIR)/mod_numfilter.o      : mod_numfilter.f90      mod_adm.o mod_grd.o mod_time.o mod_gmtr.o mod_gm_statistics.o mod_oprt.o mod_oprt3d.o mod_comm.o mod_src.o mod_vmtr.o mod_runconf.o mod_bsstate.o
$(BUILD_DIR)/mod_src_tracer.o     : mod_src_tracer.f90     mod_adm.o mod_grd.o mod_comm.o mod_gmtr.o mod_vmtr.o
$(BUILD_DIR)/mod_src.o            : mod_src.f90            mod_adm.o mod_grd.o mod_vmtr.o mod_oprt.o mod_runconf.o mod_gtl.o
$(BUILD_DIR)/mod_vi.o             : mod_vi.f90             mod_adm.o mod_grd.o mod_comm.o mod_time.o mod_oprt.o mod_vmtr.o mod_runconf.o mod_bndcnd.o mod_numfilter.o mod_src.o mod_cnvvar.o mod_bsstate.o
$(BUILD_DIR)/mod_dynamics.o       : mod_dynamics.f90       mod_adm.o mod_grd.o mod_time.o mod_vmtr.o mod_runconf.o mod_bsstate.o mod_bndcnd.o mod_prgvar.o mod_thrmdyn.o mod_numfilter.o mod_vi.o mod_src.o mod_src_tracer.o mod_nudge.o mod_oprt.o
$(BUILD_DIR)/mod_forcing_driver.o : mod_forcing_driver.f90 mod_adm.o mod_grd.o mod_gmtr.o mod_vmtr.o mod_time.o mod_runconf.o mod_cnvvar.o mod_prgvar.o mod_thrmdyn.o mod_gtl.o mod_bsstate.o mod_bndcnd.o mod_af_heldsuarez.o mod_history.o mod_af_trcadv.o
$(BUILD_DIR)/mod_gm_driver.o      : mod_gm_driver.f90      mod_adm.o mod_fio.o mod_comm.o mod_grd.o mod_gmtr.o mod_oprt.o mod_vmtr.o mod_time.o mod_extdata.o mod_runconf.o mod_prgvar.o mod_dynamics.o mod_forcing_driver.o mod_history.o mod_history_vars.o mod_embudget.o mod_bsstate.o mod_numfilter.o mod_vi.o
$(BUILD_DIR)/mod_af_heldsuarez.o  : mod_af_heldsuarez.f90  mod_adm.o
$(BUILD_DIR)/mod_af_trcadv.o      : mod_af_trcadv.f90      mod_adm.o
$(BUILD_DIR)/mod_nudge.o          : mod_nudge.f90          mod_adm.o mod_time.o mod_grd.o mod_vmtr.o mod_history.o mod_gmtr.o mod_oprt.o mod_runconf.o mod_gtl.o mod_comm.o mod_extdata.o
$(BUILD_DIR)/mod_ideal_init.o     : mod_ideal_init.f90     mod_adm.o mod_runconf.o mod_vmtr.o mod_fio.o mod_time.o mod_cnvvar.o mod_vintrpl.o mod_comm.o mod_thrmdyn.o mod_vmtr.o mod_src.o dcmip_initial_conditions_test_1_2_3_v5.o
$(BUILD_DIR)/mod_history_vars.o   : mod_history_vars.f90   mod_adm.o mod_history.o mod_runconf.o mod_gmtr.o mod_gm_statistics.o mod_vmtr.o mod_prgvar.o mod_bsstate.o mod_cnvvar.o mod_thrmdyn.o mod_bndcnd.o mod_grd.o
$(BUILD_DIR)/mod_mkgrd.o          : mod_mkgrd.f90          mod_adm.o mod_grd.o mod_comm.o mod_gm_statistics.o mod_gmtr.o
$(BUILD_DIR)/mod_mnginfo_light.o  : mod_mnginfo_light.f90  mod_misc.o
$(BUILD_DIR)/mod_option.o         : mod_option.f90
$(BUILD_DIR)/mod_netcdf.o         : mod_netcdf.f90

$(BUILD_DIR)/dcmip_initial_conditions_test_1_2_3_v5.o : dcmip_initial_conditions_test_1_2_3_v5.f90 mod_adm.o



allclean: distclean
	$(MAKE) -C $(DCUTILSDIR)      allclean
	$(MAKE) -C $(GTOOLDIR)        allclean
	$(MAKE) -C $(SCALELIBDIR)/src allclean
	rm -f $(BINDIR)/$(BINNAME)
	rm -f $(BINDIR)/$(PPNAME1)

distclean: clean
	rm -f $(BUILD_DIR)/$(BINNAME)
	rm -f $(BUILD_DIR)/$(PPNAME1)

clean:
	rm -f $(BUILD_DIR)/*.*



.SUFFIXES:
.SUFFIXES: .o .f90 .c .mod

$(BUILD_DIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(DCUTILSDIR) -I$(GTOOLDIR) -I$(MODDIR) -I$(SCALELIBDIR)/include -I$(SCALEGMDIR)/include $(PAPI_INCLUDE) $(NETCDF_INCLUDE) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR)
$(BUILD_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $< $(MODDIROPT) $(BUILD_DIR) -I$(SCALEGMDIR)/include

%.mod: %.f90
	make $(patsubst %.f90,%.o,$<)


.PHONY: clean dcutils gtool scalelib



# parameters for fixed index
glevel  =  5
rlevel  =  0
nmpi    = 10
zlayer  = 38
diamond = 10

fixedindex:
	@cd ./include; ruby ./setindex.rb $(glevel) $(rlevel) $(nmpi) $(zlayer) $(diamond)
