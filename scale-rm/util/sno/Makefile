################################################################################
#
# Makefile
#
################################################################################

TOPDIR     = $(abspath ../../..)
SYSDEP_DIR = ../../../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude


BINNAME = sno

BINS = $(BINDIR)/$(BINNAME)$(POSTFIX)
LIBS = $(LIBDIR)/libscale$(POSTFIX).a

ifneq ($(SCALE_IGNORE_SRCDEP),T)
  DEPENDLIB = $(LIBS)
endif

VPATH = $(BUILD_DIR): \
        base:         \
        plugin:

VERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(VERSION),)
  VERSION  = $(shell cat VERSION)
else
  VERSION := $(VERSION)
endif

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
	$(MAKE) build
	$(MAKE) install

build:
	$(MAKE) info
	$(MAKE) makelib
	@echo;echo "Entering SNO...";echo "Current version is " $(VERSION)
	mkdir -p $(BUILD_DIR)
	$(MAKE) makebin
	@echo "Complete making SNO."

install:
	mkdir -p $(BINDIR)
	$(MAKE) $(BINS)

info:
	@$(MAKE) -C $(SCALELIBDIR)/src --no-print-directory info
	@echo
	@echo "SNO"
	@$(MAKE) --no-print-directory conflog

makelib:
	$(MAKE) -C $(SCALELIBDIR)/src

makebin: $(BUILD_DIR)/$(BINNAME)$(POSTFIX)

$(BUILD_DIR)/$(BINNAME)$(POSTFIX) : $(BUILD_DIR)/prg_$(BINNAME).o  $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS) $(CONTRIB_LIBS)

$(BINDIR)/$(BINNAME)$(POSTFIX)  : $(BUILD_DIR)/$(BINNAME)$(POSTFIX)
	install $< $@


ifneq ($(SCALE_IGNORE_SRCDEP),T)
  $(BUILD_DIR)/prg_$(BINNAME).o : $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
endif


allclean: distclean
	$(MAKE) -C $(SCALELIBDIR)/src allclean
	rm -f $(BINDIR)/$(BINNAME)*
	rm -rf $(TOPDIR)/bin

distclean: clean
	rm -f $(BINDIR)/$(BINNAME)$(POSTFIX)
	rm -rf ./.libs*

clean:
	rm -rf $(BUILD_DIR)
	rm -f *.o *.mod *.lst *~

dependency:
	./makedepend .


.SUFFIXES:
.SUFFIXES: .o .F90 .mod

%.mod: %.F90
	make $(patsubst %.F90,%.o,$<)

$(BUILD_DIR)/%.o : %.F90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(MODDIR) $(CONTRIB_INCLUDE) -I$(SCALELIBDIR)/include -I$(SCALERMDIR)/include $(SCALE_PAPI_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

.PHONY : clean distclean allclean scalelib depend dependency

include $(TOPDIR)/utils/make/Make_environments

ifneq ($(SCALE_IGNORE_SRCDEP),T)
-include depend
endif
