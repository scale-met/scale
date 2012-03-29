################################################################################
#
# Makefile for main program
#
################################################################################

include ../../sysdep/Makedef.$(SCALE_SYS)

VPATH = ../../include:../../common:../../io:../../grid:../../atmos:../../ocean

PRJ1 = init_coldbubble
PRJ2 = init_warmbubble
PRJ3 = init_supercell
PRJ4 = init_turbdyn
PRJ5 = init_lambwave
PRJ6 = init_advection
PRJ7 = init_diffusion
PRJ8 = init_cwc

PRJS = $(PRJ1) $(PRJ2) $(PRJ3) $(PRJ4) $(PRJ5) $(PRJ6) $(PRJ7) $(PRJ8)

MODS =	\
	mod_stdio.o	\
	mod_process.o	\
	mod_const.o	\
	mod_time.o	\
	mod_grid_cartesian_plane.o	\
	mod_fileio.o	\
	mod_fileio_h.o	\
	fio.o	\
	fiof.o	\
	mod_comm.o	\
	mod_atmos_vars.o	\
	mod_atmos_sub_hydrostatic.o	\
	mod_atmos_sub_refstate.o


all: makedir makebin

makedir:
	mkdir -p ../../../bin/$(SCALE_SYS)

makebin: $(PRJS)
	install $(PRJS) ../../../bin/$(SCALE_SYS)

$(PRJ1): $(PRJ1).o
	$(LD) $(LFLAGS) -o $(PRJ1) $(PRJ1).o ${MODS}

$(PRJ2): $(PRJ2).o
	$(LD) $(LFLAGS) -o $(PRJ2) $(PRJ2).o ${MODS}

$(PRJ3): $(PRJ3).o
	$(LD) $(LFLAGS) -o $(PRJ3) $(PRJ3).o ${MODS}

$(PRJ4): $(PRJ4).o
	$(LD) $(LFLAGS) -o $(PRJ4) $(PRJ4).o ${MODS}

$(PRJ5): $(PRJ5).o
	$(LD) $(LFLAGS) -o $(PRJ5) $(PRJ5).o ${MODS}

$(PRJ6): $(PRJ6).o
	$(LD) $(LFLAGS) -o $(PRJ6) $(PRJ6).o ${MODS}

$(PRJ7): $(PRJ7).o
	$(LD) $(LFLAGS) -o $(PRJ7) $(PRJ7).o ${MODS}

$(PRJ8): $(PRJ8).o
	$(LD) $(LFLAGS) -o $(PRJ8) $(PRJ8).o ${MODS}

$(PRJ1).o: $(PRJ1).f90 mod_stdio.o ${MODS}
$(PRJ2).o: $(PRJ2).f90 mod_stdio.o ${MODS}
$(PRJ3).o: $(PRJ3).f90 mod_stdio.o ${MODS}
$(PRJ4).o: $(PRJ4).f90 mod_stdio.o ${MODS}
$(PRJ5).o: $(PRJ5).f90 mod_stdio.o ${MODS}
$(PRJ6).o: $(PRJ6).f90 mod_stdio.o ${MODS}
$(PRJ7).o: $(PRJ7).f90 mod_stdio.o ${MODS}
$(PRJ8).o: $(PRJ8).f90 mod_stdio.o ${MODS}

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
mod_atmos_sub_hydrostatic.o	: mod_atmos_sub_hydrostatic.f90 mod_stdio.o mod_const.o mod_grid_cartesian_plane.o
mod_atmos_sub_refstate.o	: mod_atmos_sub_refstate.f90 mod_stdio.o mod_process.o mod_const.o mod_grid_cartesian_plane.o mod_fileio.o mod_fileio_h.o

clean:
	rm -f $(PRJS) *.mod *.o *.lst

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.f90.o:
	$(FC) $(FFLAGS) -I../../include -o $@ -c $<
.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

%.o: %.mod
