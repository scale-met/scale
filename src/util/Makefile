################################################################################
#
# Makefile for frontend utility
#
################################################################################

include ../sysdep/Makedef.$(SCALE_SYS)

VPATH = ../common:../io:../grid:../atmos:../ocean

PRJ1 = spddump
PRJ2 = spd2bin
PRJ3 = check_diff
PRJ4 = bin2txt
PRJ5 = spddiff

PRJS = $(PRJ1) $(PRJ2) $(PRJ3) $(PRJ4) $(PRJ5)

MODS =	\
	mod_stdio.o	\
	mod_process_nompi.o	\
	mod_const.o	\
	mod_time.o	\
	fio.o	\
	fiof.o	\
	mod_fileio_h.o	\
	mod_option.o

all: makedir makebin

makedir:
	mkdir -p ../../sbin/$(SCALE_SYS)

makebin: $(PRJS)
	install $(PRJS) ../../sbin/$(SCALE_SYS)

$(PRJ1): prg_$(PRJ1).o
	$(FRONT_LD) $(FRONT_LFLAGS) -o $(PRJ1) prg_$(PRJ1).o mod_fileio_h.o fiof.o fio.o

$(PRJ2): prg_$(PRJ2).o
	$(FRONT_LD) $(FRONT_LFLAGS) -o $(PRJ2) prg_$(PRJ2).o mod_stdio.o mod_const.o mod_time.o mod_process_nompi.o mod_fileio_h.o mod_option.o fiof.o fio.o

$(PRJ3): prg_$(PRJ3).o
	$(FRONT_LD) $(FRONT_LFLAGS) -o $(PRJ3) prg_$(PRJ3).o

$(PRJ4): prg_$(PRJ4).o
	$(FRONT_LD) $(FRONT_LFLAGS) -o $(PRJ4) prg_$(PRJ4).o

$(PRJ5): prg_$(PRJ5).o
	$(FRONT_LD) $(FRONT_LFLAGS) -o $(PRJ5) prg_$(PRJ5).o mod_stdio.o mod_const.o mod_time.o mod_process_nompi.o mod_fileio_h.o mod_option.o fiof.o fio.o


mod_stdio.o         : mod_stdio.f90
mod_process_nompi.o : mod_process_nompi.f90 mod_stdio.o
mod_const.o         : mod_const.f90 mod_stdio.o mod_process_nompi.o
mod_time.o          : mod_time.f90 mod_stdio.o mod_process_nompi.o
mod_fileio_h.o	     : mod_fileio_h.f90
fio.o               : fio.c fio.h fio_def.h
fiof.o              : fiof.c fiof.h fio.h fio_def.h
mod_option.o        : mod_option.f90

prg_$(PRJ1).o	: prg_$(PRJ1).f90 mod_fileio_h.o fiof.o fio.o
prg_$(PRJ2).o	: prg_$(PRJ2).f90 mod_stdio.o mod_const.o mod_time.o mod_fileio_h.o mod_option.o fiof.o fio.o
prg_$(PRJ3).o	: prg_$(PRJ3).f90
prg_$(PRJ4).o	: prg_$(PRJ4).f90
prg_$(PRJ5).o	: prg_$(PRJ5).f90 mod_stdio.o mod_const.o mod_time.o mod_fileio_h.o mod_option.o


.SUFFIXES: .o .f90 .f .c
.f90.o:
	$(FRONT_F90) -DNO_MPIIO $(FRONT_FFLAGS) -o $@ -c $<

.c.o:
	$(FRONT_CC) -DNO_MPIIO $(FRONT_CFLAGS) -o $@ -c $<

%.o: %.mod


.PHONY: clean

clean:
	rm -f $(PRJS) *.mod *.o *.lst
