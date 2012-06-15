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


#default index, tracer, & precision set
index     = 100x10x10
tracer    = dry
precision = real8

FORMAT = netcdf

postfix = $(FORMAT)

PROJ_EXE = netcdf

FORMAT_OBJ = gtool_$(FORMAT).o

#CFLAGS_LOCAL = -DHDF5_CONFIG

all: makedir makeinclude makebin scale3_init

makedir:
	mkdir -p ../../../bin/$(SCALE_SYS)

makeinclude:
	cp ../../include/inc_precision_$(precision).f90 ../../include/inc_precision.h
	cp ../../include/inc_tracer_$(tracer).f90 ../../include/inc_tracer.h
	cp ../../include/inc_index_$(index).f90   ../../include/inc_index.h

makebin: $(PROJ_EXE)
	install $(PROJ_EXE) ../../../bin/$(SCALE_SYS)/$(PROJ_EXE)_$(index)_$(tracer)_$(precision)_$(postfix)

$(PROJ_EXE): $(PROJ_EXE).o gtool_history.o gtool_file_f.o $(FORMAT_OBJ)
	$(LD) $(LFLAGS) -o $(PROJ_EXE) $(PROJ_EXE).o $(MODS) gtool_history.o gtool_file.o gtool_file_f.o $(FORMAT_OBJ) -L/ap/netcdf/4.1.2/lib -L/ap/hdf5/1.8.7/lib -lnetcdf -lhdf5

scale3_init: scale3_init.o mod_mkinit.o mod_atmos_sub_hydrostatic.o mod_random.o
	$(LD) $(LFLAGS) -o $@ $^ $(MODS) gtool_history.o gtool_file.o gtool_file_f.o $(FORMAT_OBJ) -L/ap/netcdf/4.1.2/lib -L/ap/hdf5/1.8.7/lib -lnetcdf -lhdf5

mod_random.o: ../../common/mod_random.f90
mod_atmos_sub_hydrostatic.o: ../../atmos/mod_atmos_sub_hydrostatic.f90
mod_mkinit.o: ../../preprocess/mod_mkinit.f90 mod_random.o mod_atmos_sub_hydrostatic.o
scale3_init.o: ${PROJ_EXE} scale3_init.f90 mod_mkinit.o mod_atmos_sub_hydrostatic.o mod_random.o $(FORMAT_OBJ) gtool_file.o

gtool_file.o: gtool_file.f90 gtool_file_h.o
gtool_history.o: gtool_history.f90 gtool_file.o
gtool_file_f.o: gtool_file_f.c

mod_atmos_vars.o: gtool_file.o
mod_atmos_vars_sf.o: gtool_file.o
mod_atmos_sub_boundary.o: gtool_file.o
mod_atmos_sub_refstate.o: gtool_file.o
mod_getmetrics.o: gtool_file.o
mod_grid_cartesian_plane.o: gtool_file.o
mod_ocean_vars.o: gtool_file.o

mod_history.o: gtool_history.o


.SUFFIXES: .o .f90 .f .c
.f90.o:
	$(FC) $(FFLAGS) $(CFLAGS_LOCAL) -I../../include -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) $(CFLAGS_LOCAL) -I/ap/netcdf/4.1.2/include -o $@ -c $<

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
