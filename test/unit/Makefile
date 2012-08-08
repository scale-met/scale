TOPDIR = ../..
BINDIR = $(TOPDIR)/bin
LIBDIR = $(TOPDIR)/lib
SRCDIR = $(TOPDIR)/src
GTOOLDIR = $(SRCDIR)/gtool

include $(SRCDIR)/sysdep/Makedef.$(SCALE_SYS)

VPATH = ../..:$(SRCDIR)/common:$(SRCDIR)/io:$(SRCDIR)/grid:$(SRCDIR)/atmos:$(SRCDIR)/ocean:$(SRCDIR)/preprocess

PROJ_EXE = unit

#default index & tracer set
index   = 10x10x2
tracer  = kessler
precision = real8
postfix = 

LIBNAME = libscale_$(index)_$(tracer)_$(precision)_$(postfix).a

all: makeinclude makebin test

makeinclude:
	cp $(SRCDIR)/include/inc_tracer_$(tracer).f90 $(SRCDIR)/include/inc_tracer.h
	cp $(SRCDIR)/include/inc_index_$(index).f90   $(SRCDIR)/include/inc_index.h
	cp $(SRCDIR)/include/inc_precision_$(precision).f90 $(SRCDIR)/include/inc_precision.h

makebin: $(PROJ_EXE)

test: makebin
	./$(PROJ_EXE) dummy.cnf

$(PROJ_EXE): $(PROJ_EXE).o test_atmos_phy_tb_smg.o test_atmos_dyn_fent_fct.o dc_test.o dc_types.o $(LIBDIR)/$(LIBNAME) $(LIBDIR)/libgtool.a
	$(LD) $(LFLAGS) -o $@ $^ $(NETCDF_LIBS)

$(LIBDIR)/libgtool.a:
	cd $(GTOOLDIR); make

$(LIBDIR)/$(LIBNAME):
	cd $(SRCDIR); make distclean; make index=$(index) tracer=$(tracer) precision=$(precision)

$(PROJ_EXE).o: $(PROJ_EXE).f90 test_atmos_phy_tb_smg.o test_atmos_dyn_fent_fct.o $(LIBDIR)/libgtool.a $(LIBDIR)/$(LIBNAME)

test_atmos_phy_tb_smg.o: test_atmos_phy_tb_smg.f90 dc_test.o $(LIBDIR)/libgtool.a $(LIBDIR)/$(LIBNAME)
test_atmos_dyn_fent_fct.o: test_atmos_dyn_fent_fct.f90 dc_test.o $(LIBDIR)/libgtool.a $(LIBDIR)/$(LIBNAME)

dc_test.o: dc_types.o

.SUFFIXES: .o .f90 .f .c
.f90.o:
	$(FC) -DDEBUG $(FFLAGS) -O0 -g -I$(SRCDIR)/include -I$(SRCDIR) -I$(GTOOLDIR) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) -O0 -g -o $@ -c $<

%.o: %.mod


.PHONY: clean distclean

clean:
	rm -f *.mod *.o *~ *.lst

distclean: clean
	rm -f $(PROJ_EXE)
