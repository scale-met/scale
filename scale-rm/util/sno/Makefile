################################################################################
#
# Makefile
#
################################################################################

TOPDIR     = $(abspath ../../..)
BUILD_DIR  = ./.libs
SYSDEP_DIR = ../../../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

BINNAME = sno

LIBS = $(LIBDIR)/libscale*.a

VPATH = plugin:

OBJS = mod_sno_h.o     \
       mod_sno.o       \
       mod_sno_map.o   \
       mod_sno_axis.o  \
       mod_sno_grads.o \
       mod_sno_comm.o  \
       mod_sno_vars.o  \
       \
       mod_snoplugin_timeave.o \
       mod_snoplugin_hgridope.o \
       mod_snoplugin_vgridope.o

all:
	$(MAKE) envlog
	$(MAKE) makedir
	$(MAKE) makebin
	@echo "Complete making."

makedir:
	mkdir -p $(BINDIR)
	mkdir -p $(BUILD_DIR)

makebin: $(BINDIR)/$(BINNAME)

$(BINDIR)/$(BINNAME) : $(BUILD_DIR)/prg_$(BINNAME).o  $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(LIBS) $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS) $(CONTRIB_LIBS)

$(BUILD_DIR)/prg_$(BINNAME).o : prg_$(BINNAME).F90 $(patsubst %,$(BUILD_DIR)/%,$(OBJS))

$(BUILD_DIR)/mod_sno_h.o     : mod_sno_h.F90
$(BUILD_DIR)/mod_sno.o       : mod_sno.F90       $(BUILD_DIR)/mod_sno_h.o
$(BUILD_DIR)/mod_sno_map.o   : mod_sno_map.F90   $(BUILD_DIR)/mod_sno_h.o
$(BUILD_DIR)/mod_sno_axis.o  : mod_sno_axis.F90  $(BUILD_DIR)/mod_sno_h.o $(BUILD_DIR)/mod_sno.o
$(BUILD_DIR)/mod_sno_grads.o : mod_sno_grads.F90 $(BUILD_DIR)/mod_sno_h.o
$(BUILD_DIR)/mod_sno_comm.o  : mod_sno_comm.F90  $(BUILD_DIR)/mod_sno_h.o $(BUILD_DIR)/mod_sno.o $(BUILD_DIR)/mod_sno_axis.o
$(BUILD_DIR)/mod_sno_vars.o  : mod_sno_vars.F90  $(BUILD_DIR)/mod_sno_h.o $(BUILD_DIR)/mod_sno.o $(BUILD_DIR)/mod_sno_axis.o $(BUILD_DIR)/mod_sno_grads.o $(BUILD_DIR)/mod_sno_comm.o

$(BUILD_DIR)/mod_snoplugin_timeave.o  : mod_snoplugin_timeave.F90  $(BUILD_DIR)/mod_sno_h.o $(BUILD_DIR)/mod_sno_vars.o
$(BUILD_DIR)/mod_snoplugin_hgridope.o : mod_snoplugin_hgridope.F90 $(BUILD_DIR)/mod_sno_h.o $(BUILD_DIR)/mod_sno_vars.o
$(BUILD_DIR)/mod_snoplugin_vgridope.o : mod_snoplugin_vgridope.F90 $(BUILD_DIR)/mod_sno_h.o $(BUILD_DIR)/mod_sno_vars.o



allclean: distclean
	rm -f $(BINDIR)/$(BINNAME)

distclean: clean
	rm -f $(BUILD_DIR)/$(BINNAME)

clean:
	rm -f $(BUILD_DIR)/*.*



.SUFFIXES:
.SUFFIXES: .o .F90 .mod

%.mod: %.F90
	make $(patsubst %.F90,%.o,$<)

$(BUILD_DIR)/%.o : %.F90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(MODDIR) $(CONTRIB_INCLUDE) -I$(SCALELIBDIR)/include -I$(SCALERMDIR)/include $(SCALE_PAPI_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

.PHONY : clean distclean allclean

include $(TOPDIR)/utils/make/Make_environments
