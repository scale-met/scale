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
index   = 436x16x16
tracer  = ndw6
postfix = 

PROJ_EXE = HDF5_full

all: makedir makeinclude makebin

makedir:
	mkdir -p ../../../bin/$(SCALE_SYS)

makeinclude:
	cp ../../include/inc_tracer_$(tracer).f90 ../../include/inc_tracer.h
	cp ../../include/inc_index_$(index).f90   ../../include/inc_index.h

makebin: $(PROJ_EXE)
	install $(PROJ_EXE) ../../../bin/$(SCALE_SYS)/$(PROJ_EXE)_$(index)_$(tracer)_$(postfix)

$(PROJ_EXE): $(PROJ_EXE).o
	$(LD) $(LFLAGS) -o $(PROJ_EXE) $(PROJ_EXE).o $(MODS) $(HDF5_LIBS)


.SUFFIXES: .o .f90 .f .c
.f90.o:
	$(FC) -DCONFIG_HDF5 -DCONFIG_HDF5_INPUT $(HDF5_INCLUDE) $(FFLAGS) -I../../include -o $@ -c $<

.c.o:
	$(CC) -DCONFIG_HDF5 -DCONFIG_HDF5_INPUT $(HDF5_INCLUDE) $(CFLAGS) -o $@ -c $<

%.o: %.mod


.PHONY: depends freeze clean copylocal mktar

clean:
	rm -f $(PROJ_EXE) depends.info *.mod *.o *~ *.lst

freeze: clean depends mktar

mktar:
	tar -czf ../$(PROJ_EXE).tgz ./

depends:
	@../mkdepends.pl $(PROJ_EXE) $(VPATH)

-include depends.info
