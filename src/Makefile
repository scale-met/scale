################################################################################
#
# Makefile for main program
#
################################################################################

include ./sysdep/Makedef.$(SCALE_SYS)

VPATH = common:io:grid:atmos:ocean:preprocess/initial

PRJ1 = scale3
PRJ2 = init_coldbubble

PRJS = $(PRJ1) $(PRJ2)

MODS =	\
	mod_stdio.o	\
	mod_process.o	\
	mod_const.o	\
	mod_time.o	\
	mod_fileio.o	\
	mod_fileio_h.o	\
	fio.o	\
	fiof.o	\
	mod_grid_cartesian_plane.o	\
	mod_comm.o	\
	mod_atmos.o	\
	mod_atmos_vars.o	\
	mod_atmos_sub_refstate.o	\
	mod_atmos_sub_boundary.o	\
	mod_atmos_dyn_fent_pdfct.o	\
	mod_atmos_phy_sf.o	\
	mod_ocean_vars.o	\
	mod_atmos_phy_tb_smg.o	\
	mod_history.o	\
	mod_atmos_sub_hydrostatic.o	\
	mod_ocean.o

all: makedir makebin

makedir:
	mkdir -p ../bin/$(SCALE_SYS)

makebin: $(PRJS)
	install $(PRJS) ../bin/$(SCALE_SYS)

$(PRJ1): $(PRJ1).o
	$(LD) $(LFLAGS) -o $(PRJ1) $(PRJ1).o ${MODS}

$(PRJ2): $(PRJ2).o
	$(LD) $(LFLAGS) -o $(PRJ2) $(PRJ2).o ${MODS}

$(PRJ1).o: $(PRJ1).f90 ${MODS}
$(PRJ2).o: $(PRJ2).f90 ${MODS}

mod_stdio.o	: mod_stdio.f90
mod_process.o	: mod_process.f90 mod_stdio.o
mod_const.o	: mod_const.f90 mod_stdio.o mod_process.o
mod_time.o	: mod_time.f90 mod_stdio.o mod_process.o
mod_fileio_h.o	: mod_fileio_h.f90
fio.o	: fio.c fio.h fio_def.h
fiof.o	: fiof.c fiof.h fio.h fio_def.h
mod_fileio.o	: mod_fileio.f90 mod_fileio_h.o mod_stdio.o mod_process.o mod_time.o mod_const.o fio.o fiof.o
mod_grid_cartesian_plane.o	: mod_grid_cartesian_plane.f90 mod_stdio.o mod_process.o mod_fileio.o mod_fileio_h.o
mod_comm.o	: mod_comm.f90 mod_stdio.o mod_process.o mod_grid_cartesian_plane.o mod_time.o mod_const.o
mod_atmos_vars.o	: mod_atmos_vars.f90 mod_stdio.o mod_fileio_h.o mod_process.o mod_const.o mod_grid_cartesian_plane.o mod_comm.o mod_fileio.o mod_time.o
mod_atmos_sub_refstate.o	: mod_atmos_sub_refstate.f90 mod_stdio.o mod_process.o mod_const.o mod_grid_cartesian_plane.o mod_fileio.o mod_fileio_h.o
mod_atmos_sub_boundary.o	: mod_atmos_sub_boundary.f90 mod_stdio.o mod_fileio_h.o mod_process.o mod_const.o mod_grid_cartesian_plane.o mod_comm.o mod_fileio.o
mod_atmos_dyn_fent_pdfct.o	: mod_atmos_dyn_fent_pdfct.f90 mod_stdio.o mod_process.o mod_const.o mod_grid_cartesian_plane.o mod_time.o mod_comm.o mod_atmos_vars.o mod_atmos_sub_refstate.o mod_atmos_sub_boundary.o
mod_ocean_vars.o	: mod_ocean_vars.f90 mod_stdio.o mod_fileio_h.o mod_process.o mod_const.o mod_grid_cartesian_plane.o mod_comm.o mod_fileio.o mod_time.o
mod_atmos_phy_sf.o	: mod_atmos_phy_sf.f90 mod_stdio.o mod_process.o mod_const.o mod_grid_cartesian_plane.o mod_atmos_vars.o mod_ocean_vars.o
mod_atmos_phy_tb_smg.o	: mod_atmos_phy_tb_smg.f90 mod_const.o mod_grid_cartesian_plane.o mod_atmos_vars.o
mod_history.o	: mod_history.f90 mod_stdio.o mod_fileio_h.o mod_process.o mod_time.o mod_grid_cartesian_plane.o mod_fileio.o
mod_atmos.o	: mod_atmos.f90 mod_stdio.o mod_atmos_vars.o mod_atmos_sub_refstate.o mod_atmos_sub_boundary.o mod_atmos_dyn_fent_pdfct.o mod_time.o mod_grid_cartesian_plane.o mod_atmos_phy_sf.o mod_atmos_phy_tb_smg.o mod_history.o
mod_ocean.o	: mod_ocean.f90 mod_stdio.o mod_ocean_vars.o mod_grid_cartesian_plane.o
mod_atmos_sub_hydrostatic.o	: mod_atmos_sub_hydrostatic.f90 mod_stdio.o mod_const.o mod_grid_cartesian_plane.o

clean:
	rm -f $(PRJS) *.mod *.o *.lst

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.f90.o:
	$(FC) $(FFLAGS) -o $@ -c $<
.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

%.o: %.mod
