################################################################################
#
# Makefile for main program
#
################################################################################

include ../../sysdep/Makedef.$(SCALE_SYS)

VPATH = ../../common:../../io:../../grid:../../atmos:../../ocean

PRJ1 = init_coldbubble
PRJ2 = init_warmbubble

PRJS = $(PRJ1) $(PRJ2)

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
	mod_atmos_vars.o

all: makedir makebin

makedir:
	mkdir -p ../../../bin/$(SCALE_SYS)

makebin: $(PRJS)
	install $(PRJS) ../../../bin/$(SCALE_SYS)

$(PRJ1): $(PRJ1).o
	$(LD) $(LFLAGS) -o $(PRJ1) $(PRJ1).o ${MODS}

$(PRJ2): $(PRJ2).o
	$(LD) $(LFLAGS) -o $(PRJ2) $(PRJ2).o ${MODS}

$(PRJ1).o: $(PRJ1).f90 mod_stdio.o mod_process.o mod_const.o mod_time.o mod_grid_cartesian_plane.o mod_comm.o mod_fileio.o mod_atmos_vars.o
$(PRJ2).o: $(PRJ2).f90 mod_stdio.o mod_process.o mod_const.o mod_time.o mod_grid_cartesian_plane.o mod_comm.o mod_fileio.o mod_atmos_vars.o

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

clean:
	rm -f $(PRJS) *.mod *.o

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.f90.o:
	$(FC) $(FFLAGS) -o $@ -c $<
.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
