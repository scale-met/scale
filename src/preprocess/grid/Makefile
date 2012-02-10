################################################################################
#
# Makefile for main program
#
################################################################################

include ../../sysdep/Makedef.$(SCALE_SYS)

VPATH = ../../common:../../io:../../grid:../../atmos:../../ocean

PRJ1 = makegrid

PRJS = $(PRJ1)

MODS =	\
	mod_stdio.o	\
	mod_process.o	\
	mod_const.o	\
	mod_time.o	\
	mod_fileio.o	\
	mod_fileio_h.o	\
	fio.o	\
	fiof.o	\
	mod_grid_cartesian_plane.o

all: makedir makebin

makedir:
	mkdir -p ../../../bin/$(SCALE_SYS)

makebin: $(PRJS)
	install $(PRJS) ../../../bin/$(SCALE_SYS)

$(PRJ1): $(PRJ1).o
	$(LD) $(LFLAGS) -o $(PRJ1) $(PRJ1).o ${MODS}

$(PRJ1).o: $(PRJ1).f90 mod_stdio.o mod_process.o mod_const.o mod_time.o mod_fileio.o mod_grid_cartesian_plane.o

mod_stdio.o	: mod_stdio.f90
mod_process.o	: mod_process.f90 mod_stdio.o
mod_const.o	: mod_const.f90 mod_stdio.o mod_process.o
mod_time.o	: mod_time.f90 mod_stdio.o mod_process.o
mod_fileio_h.o	: mod_fileio_h.f90
fio.o	: fio.c fio.h fio_def.h
fiof.o	: fiof.c fiof.h fio.h fio_def.h
mod_fileio.o	: mod_fileio.f90 mod_fileio_h.o mod_stdio.o mod_process.o mod_time.o mod_const.o fio.o fiof.o
mod_grid_cartesian_plane.o	: mod_grid_cartesian_plane.f90 mod_stdio.o mod_process.o mod_fileio.o mod_fileio_h.o

clean:
	rm -f $(PRJS) *.mod *.o *.lst

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.f90.o:
	$(FC) $(FFLAGS) -o $@ -c $<
.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
