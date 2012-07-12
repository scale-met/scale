################################################################################
#
# Makefile for each test project
# %> make clean       (cleanup)
# %> make depends     (generate depends.info and include)
# %> make
#
# if you have the problem about compilation order,
# please try the following commands;
#
# %> make -n       (display all procedure for check)
# or
# %> make -B       (ignore the timestamp)
#
################################################################################
include ../../sysdep/Makedef.$(SCALE_SYS)

VPATH = ../../common:../../io:../../grid:../../atmos:../../ocean:../../preprocess

#default index & tracer set
index   = 60x45x45
#index   = 60x90x90
tracer  = dry
precision = real8
postfix = 

PROJ_EXE = smagorinsky

all: makedir makeinclude makebin

makedir:
	mkdir -p ../../../bin/$(SCALE_SYS)

makeinclude:
	cp ../../include/inc_tracer_$(tracer).f90 ../../include/inc_tracer.h
	cp ../../include/inc_index_$(index).f90   ../../include/inc_index.h
	cp ../../include/inc_precision_$(precision).f90 ../../include/inc_precision.h

makebin: $(PROJ_EXE)
	install $(PROJ_EXE) ../../../bin/$(SCALE_SYS)/$(PROJ_EXE)_$(index)_$(tracer)_$(postfix)

$(PROJ_EXE): $(PROJ_EXE).o
#	$(LD) -g -O0 -fpe0 -traceback -o $(PROJ_EXE) $(PROJ_EXE).o $(MODS)
	$(LD) $(LFLAGS) -o $(PROJ_EXE) $(PROJ_EXE).o $(MODS)

scale3_init: scale3_init.o
	$(LD) $(LFLAGS) -o $@ $< mod_random.o mod_mkinit.o mod_atmos_sub_hydrostatic.o $(MODS)

scale3_init.o: scale3_init.f90 mod_random.o mod_mkinit.o
	$(FC) $(FFLAGS) -I../../include -o $@ -c $<

mod_mkinit.o: mod_atmos_sub_hydrostatic.o


.SUFFIXES: .o .f90 .f .c
.f90.o:
#	$(FC) -fpe0 -g -O0 -fpp3 -traceback -I../../include -o $@ -c $<
	$(FC)  $(FFLAGS) -I../../include -o $@ -c $<

.c.o:
	$(CC)  $(CFLAGS) -o $@ -c $<

%.o: %.mod


.PHONY: depends freeze clean copylocal mktar

clean:
	rm -f $(PROJ_EXE) depends.info *.mod *.o *~ *.lst

freeze: clean depends copylocal mktar

copylocal:
	@for i in $(MODSRC) ; do \
		cp $$i . ; \
	done
	cp ../../sysdep/Makedef.$(SCALE_SYS) .

mktar:
	tar -czf ../$(PROJ_EXE).tgz ./

depends:
	@../mkdepends.pl $(PROJ_EXE) $(VPATH)

-include depends.info
