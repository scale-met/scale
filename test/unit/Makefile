TOPDIR = ../..
LIBDIR = $(TOPDIR)/lib
SRCDIR = $(TOPDIR)/src
INCDIR = $(SRCDIR)/include

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)

VPATH = $(SRCDIR)/common:$(SRCDIR)/io:$(SRCDIR)/grid:$(SRCDIR)/atmos:$(SRCDIR)/ocean:$(SRCDIR)/preprocess

PROJ_EXE = unit

#default index & tracer set
index   = 10x10x2
tracer  = kessler
sfcflx  = const
precision = real8
postfix = debug

LIBNAME = libscale_$(index)_$(tracer)_$(sfcflx)_$(precision)_$(postfix).a

MODS = \
	test_atmos_phy_tb_smg.o \
	test_atmos_dyn_fent_fct.o \
	dc_test.o \
	dc_types.o

all: makelib makeinclude makebin

makelib:
	cd $(SRCDIR) && $(MAKE) index=$(index) tracer=$(tracer) sfcflx=$(sfcflx) precision=$(precision) postfix=$(postfix) FFLAGS="$(FFLAGS_DEBUG)"
#	cd $(SRCDIR) && $(MAKE) clean && $(MAKE) index=$(index) tracer=$(tracer) precision=$(precision) postfix=$(postfix) FFLAGS="$(FFLAGS_DEBUG)"

makeinclude:
	cp $(INCDIR)/inc_tracer_$(tracer).f90       $(INCDIR)/inc_tracer.h
	cp $(INCDIR)/inc_index_$(index).f90         $(INCDIR)/inc_index.h
	cp $(INCDIR)/inc_precision_$(precision).f90 $(INCDIR)/inc_precision.h

makebin: $(PROJ_EXE)

$(PROJ_EXE): $(PROJ_EXE).o $(MODS) $(LIBDIR)/$(LIBNAME)
	$(LD) $(LFLAGS) -o $(PROJ_EXE) $^

$(PROJ_EXE).o: $(MODS)

test_atmos_phy_tb_smg.o: test_atmos_phy_tb_smg.f90 dc_test.o

test_atmos_dyn_fent_fct.o: test_atmos_dyn_fent_fct.f90 dc_test.o

dc_test.o: dc_test.f90 dc_types.o

.SUFFIXES: .o .f90 .f .c
.f90.o:
	$(FC) $(FFLAGS_DEBUG) -I$(INCDIR) -I$(SRCDIR) -o $@ -c $<
#	$(FC) -DDEBUG -O0 -g -fpp3 -fpe0 -traceback -I$(SRCDIR) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) -O0 -g -o $@ -c $<

%.o: %.mod


.PHONY: clean distclean

clean:
	rm -f *.mod *.o *~ *.lst

distclean: clean
	rm -f $(PROJ_EXE)
