TOPDIR = ../..

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/src/Mkinclude

VPATH = $(SRCDIR)/common:$(SRCDIR)/io:$(SRCDIR)/grid:$(SRCDIR)/atmos:$(SRCDIR)/ocean:$(SRCDIR)/preprocess

PROJ_EXE = unit

#default index & tracer set
index   = 10x10x2
tracer  = kessler
radiation = mstrnX
sfcflx  = const
precision = real8
postfix = debug

MODS = \
	test_atmos_phy_tb_smg.o \
	test_atmos_dyn_fent_fct.o

all: makelib makeinclude makebin run

makelib:
	cd $(SRCDIR) && $(MAKE) index=$(index) tracer=$(tracer) radiation=$(radiation) sfcflx=$(sfcflx) precision=$(precision) postfix=$(postfix) FFLAGS="$(FFLAGS_DEBUG)"

makeinclude:
	cp $(INCDIR)/inc_tracer_$(tracer).f90       $(INCDIR)/inc_tracer.h
	cp $(INCDIR)/inc_index_$(index).f90         $(INCDIR)/inc_index.h
	cp $(INCDIR)/inc_precision_$(precision).f90 $(INCDIR)/inc_precision.h

makebin: $(PROJ_EXE)

run:
	./$(PROJ_EXE) dummy.conf

$(PROJ_EXE): $(PROJ_EXE).o $(MODS) $(LIBDIR)/$(LIBNAME) $(LIBDIR)/libgtool.a $(LIBDIR)/libdcutils.a
	$(LD) $(LFLAGS) -o $(PROJ_EXE) $^ $(NETCDF_LIBS)

$(PROJ_EXE).o: $(MODS)

test_atmos_phy_tb_smg.o: test_atmos_phy_tb_smg.f90

test_atmos_dyn_fent_fct.o: test_atmos_dyn_fent_fct.f90

.SUFFIXES: .o .f90 .f
.f90.o:
	$(FC) $(FFLAGS_DEBUG) -I$(INCDIR) -I$(SRCDIR) -I$(DCUTILSDIR) -o $@ -c $<

%.o: %.mod

.PHONY: clean distclean

clean:
	rm -f *.mod *.o *~ *.lst

distclean: clean
	$(MAKE) -C $(SRCDIR) distclean
	rm -f $(PROJ_EXE)
