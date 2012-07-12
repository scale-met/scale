TOPDIR = ../../../../..
BINDIR = $(TOPDIR)/bin
SRCDIR = $(TOPDIR)/src

include $(SRCDIR)/sysdep/Makedef.$(SCALE_SYS)

VPATH = ../..:$(SRCDIR)/common:$(SRCDIR)/io:$(SRCDIR)/grid:$(SRCDIR)/atmos:$(SRCDIR)/ocean:$(SRCDIR)/preprocess

PROJ_EXE = unit

#default index & tracer set
index   = 10x10x2
tracer  = ndw6
precision = real8
postfix = 

all: makedir makeinclude makebin

makedir:
	mkdir -p $(BINDIR)/$(SCALE_SYS)

makeinclude:
	cp $(SRCDIR)/include/inc_tracer_$(tracer).f90 $(SRCDIR)/include/inc_tracer.h
	cp $(SRCDIR)/include/inc_index_$(index).f90   $(SRCDIR)/include/inc_index.h
	cp $(SRCDIR)/include/inc_precision_$(precision).f90 $(SRCDIR)/include/inc_precision.h

makebin: $(PROJ_EXE)

$(PROJ_EXE): $(PROJ_EXE).o
	$(LD) $(LFLAGS) -o $(PROJ_EXE) $(PROJ_EXE).o $(MODS) dc_test.o dc_types.o


test_atmos_phy_tb_smg.o: test_atmos_phy_tb_smg.f90 mod_atmos_phy_tb_smg.o dc_test.o

dc_test.o: dc_test.f90 dc_types.o

.SUFFIXES: .o .f90 .f .c
.f90.o:
	$(FC) -DDEBUG $(FFLAGS) -O0 -g -I$(SRCDIR)/include -o $@ -c $<
#	$(FC) -DDEBUG -O0 -g -fpp3 -fpe0 -traceback -I$(SRCDIR)/include -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) -O0 -g -o $@ -c $<

%.o: %.mod


.PHONY: depends freeze clean copylocal mktar

clean:
	rm -f *.mod *.o *~ *.lst

cleanall: clean
	rm -f $(PROJ_EXE) depends.info

freeze: clean depends copylocal mktar

copylocal:
	@for i in $(MODSRC) ; do \
		cp $$i . ; \
	done
	cp $(SRCDIR)/sysdep/Makedef.$(SCALE_SYS) .

depends:
	@../../../mkdepends.pl unit $(VPATH)

-include depends.info
