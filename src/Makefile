################################################################################
#
# Makefile for main program
#
################################################################################

include ./sysdep/Makedef.$(SCALE_SYS)

VPATH = include:common:io:grid:atmos:ocean:preprocess/initial

#default index & tracer set
index   = 336x63x63
tracer  = ndw6
postfix = 

PRJ1 = scale3
PRJ2 = init_coldbubble
PRJ3 = init_warmbubble
PRJ4 = init_turbdyn

PRJS = $(PRJ1) $(PRJ2) $(PRJ3) $(PRJ4)

MODS =	\
	mod_stdio.o	\
	mod_process.o	\
	mod_const.o	\
	mod_random.o	\
	mod_time.o	\
	mod_fileio.o	\
	mod_fileio_h.o	\
	fio.o	\
	fiof.o	\
	mod_grid_cartesian_plane.o	\
	mod_geometrics.o	\
	mod_comm.o	\
	mod_history.o	\
	mod_monitor.o	\
	mod_atmos.o	\
	mod_atmos_vars.o	\
	mod_atmos_sub_thermodyn.o	\
	mod_atmos_sub_refstate.o	\
	mod_atmos_sub_boundary.o	\
	mod_atmos_dyn_fent_pdfct.o	\
	mod_atmos_phy_sf.o	\
	mod_ocean_vars.o	\
	mod_atmos_vars_sf.o	\
	mod_atmos_phy_tb_smg.o	\
	mod_ocean.o	\
	mod_ocean_fixedsst.o \
	mod_mkexp.o	\
	mod_atmos_sub_hydrostatic.o

all: makedir makeinclude makebin

makedir:
	mkdir -p ../bin/$(SCALE_SYS)

makeinclude:
	@if [ ! -f include/inc_tracer.h ]; then \
		cp include/inc_tracer_$(tracer).f90 include/inc_tracer.h; \
	fi
	@if [ ! -f include/inc_index.h ]; then \
		cp include/inc_index_$(index).f90 include/inc_index.h; \
	fi

makebin: $(PRJS)
	install $(PRJ1) ../bin/$(SCALE_SYS)/$(PRJ1)_$(index)_$(tracer)_$(postfix)
	install $(PRJ2) ../bin/$(SCALE_SYS)/$(PRJ2)_$(index)_$(tracer)_$(postfix)
	install $(PRJ3) ../bin/$(SCALE_SYS)/$(PRJ3)_$(index)_$(tracer)_$(postfix)
	install $(PRJ4) ../bin/$(SCALE_SYS)/$(PRJ4)_$(index)_$(tracer)_$(postfix)

$(PRJ1): $(PRJ1).o
	$(LD) $(LFLAGS) -o $(PRJ1) $(PRJ1).o ${MODS}

$(PRJ2): $(PRJ2).o
	$(LD) $(LFLAGS) -o $(PRJ2) $(PRJ2).o ${MODS}

$(PRJ3): $(PRJ3).o
	$(LD) $(LFLAGS) -o $(PRJ3) $(PRJ3).o ${MODS}

$(PRJ4): $(PRJ4).o
	$(LD) $(LFLAGS) -o $(PRJ4) $(PRJ4).o ${MODS}

$(PRJ1).o: $(PRJ1).f90 $(MODS)
$(PRJ2).o: $(PRJ2).f90 $(MODS)
$(PRJ3).o: $(PRJ3).f90 $(MODS)

mod_stdio.o	: mod_stdio.f90
mod_process.o	: mod_process.f90 mod_stdio.o
mod_const.o	: mod_const.f90 mod_stdio.o mod_process.o
mod_random.o	: mod_random.f90 mod_stdio.o mod_process.o
mod_time.o	: mod_time.f90 mod_stdio.o mod_process.o
mod_fileio_h.o	: mod_fileio_h.f90
fio.o	: fio.c fio.h fio_def.h
fiof.o	: fiof.c fiof.h fio.h fio_def.h
mod_fileio.o	: mod_fileio.f90 mod_stdio.o mod_time.o mod_fileio_h.o mod_process.o mod_const.o fio.o fiof.o
mod_grid_cartesian_plane.o	: mod_grid_cartesian_plane.f90 mod_stdio.o mod_process.o mod_fileio.o mod_fileio_h.o
mod_geometrics.o	: mod_geometrics.f90 mod_stdio.o mod_process.o mod_time.o mod_fileio_h.o mod_fileio.o mod_grid_cartesian_plane.o mod_const.o
mod_comm.o	: mod_comm.f90 mod_stdio.o mod_time.o mod_process.o mod_const.o mod_geometrics.o
mod_history.o	: mod_history.f90 mod_stdio.o mod_fileio_h.o mod_process.o mod_time.o mod_fileio.o
mod_monitor.o	: mod_monitor.f90 mod_stdio.o mod_fileio_h.o mod_process.o mod_geometrics.o
mod_atmos_sub_thermodyn.o	: mod_atmos_sub_thermodyn.f90 mod_const.o mod_time.o
mod_atmos_vars.o	: mod_atmos_vars.f90 mod_stdio.o mod_fileio_h.o mod_process.o mod_const.o mod_history.o mod_comm.o mod_fileio.o mod_time.o mod_grid_cartesian_plane.o mod_atmos_sub_thermodyn.o mod_monitor.o
mod_atmos_sub_refstate.o	: mod_atmos_sub_refstate.f90 mod_stdio.o mod_process.o mod_fileio.o mod_fileio_h.o mod_const.o mod_grid_cartesian_plane.o
mod_atmos_sub_boundary.o	: mod_atmos_sub_boundary.f90 mod_stdio.o mod_fileio_h.o mod_process.o mod_const.o mod_grid_cartesian_plane.o mod_fileio.o mod_comm.o mod_time.o mod_atmos_sub_refstate.o
mod_atmos_dyn_fent_pdfct.o	: mod_atmos_dyn_fent_pdfct.f90 mod_stdio.o mod_time.o mod_process.o mod_grid_cartesian_plane.o mod_const.o mod_comm.o mod_atmos_vars.o mod_atmos_sub_refstate.o mod_atmos_sub_boundary.o
mod_ocean_vars.o	: mod_ocean_vars.f90 mod_stdio.o mod_fileio_h.o mod_process.o mod_comm.o mod_fileio.o mod_time.o
mod_atmos_vars_sf.o	: mod_atmos_vars_sf.f90 mod_stdio.o mod_fileio_h.o mod_process.o mod_comm.o mod_fileio.o mod_time.o
mod_atmos_phy_sf.o	: mod_atmos_phy_sf.f90 mod_stdio.o mod_process.o mod_grid_cartesian_plane.o mod_const.o mod_atmos_vars.o mod_ocean_vars.o mod_atmos_vars_sf.o
mod_atmos_phy_tb_smg.o	: mod_atmos_phy_tb_smg.f90 mod_stdio.o mod_const.o mod_time.o mod_comm.o mod_grid_cartesian_plane.o mod_atmos_vars.o mod_atmos_vars_sf.o
mod_atmos.o	: mod_atmos.f90 mod_stdio.o mod_time.o mod_atmos_vars.o mod_atmos_sub_refstate.o mod_atmos_sub_boundary.o mod_atmos_sub_thermodyn.o mod_atmos_dyn_fent_pdfct.o mod_atmos_phy_sf.o mod_atmos_phy_tb_smg.o
mod_ocean_fixedsst.o	: mod_ocean_fixedsst.f90 mod_stdio.o mod_process.o mod_ocean_vars.o mod_time.o
mod_ocean.o	: mod_ocean.f90 mod_stdio.o mod_time.o mod_ocean_vars.o mod_ocean_fixedsst.o
mod_mkexp.o	: mod_mkexp.f90 mod_stdio.o mod_process.o mod_const.o mod_grid_cartesian_plane.o mod_atmos_vars.o mod_atmos_sub_hydrostatic.o mod_random.o
mod_atmos_sub_hydrostatic.o	: mod_atmos_sub_hydrostatic.f90 mod_stdio.o mod_const.o mod_grid_cartesian_plane.o

clean:
	rm -f $(PRJS) include/inc_index.h include/inc_tracer.h *.mod *.o *.lst

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.f90.o:
	$(FC) $(FFLAGS) -I./include -o $@ -c $<
.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

%.o: %.mod
