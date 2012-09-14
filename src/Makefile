################################################################################
#
# Makefile for main program
#
################################################################################

TOPDIR = ..

include ../sysdep/Makedef.$(SCALE_SYS)
include ./Mkinclude

VPATH = include:common:io:grid:atmos:ocean:preprocess


#default index & tracer set
index     = 206x24x24
tracer    = kessler
precision = real8
sfcflx    = louis
radiation = mstrnX
postfix   = 

PRJ1 = scale3
PRJ2 = scale3_init

PRJS = $(PRJ1) $(PRJ2)

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
	mod_atmos_dyn_fent_fullfct_lss.o	\
	mod_atmos_phy_sf_$(sfcflx).o	\
	mod_atmos_vars_sf.o	\
	mod_atmos_phy_tb_smg.o	\
        mod_atmos_phy_mp_$(tracer).o \
	mod_atmos_sub_precipitation.o	\
        mod_atmos_phy_rd_$(radiation).o \
	mod_ocean.o	\
	mod_ocean_vars.o	\
	mod_ocean_fixedsst.o	\
	mod_mkinit.o	\
	mod_atmos_sub_hydrostatic.o	\
	mod_atmos_mp_saturation.o	\
	mod_atmos_sub_saturation.o

all: makedir makeinclude makebin

makedir:
	mkdir -p ../bin/$(SCALE_SYS)

makeinclude:
	cp include/inc_tracer_$(tracer).f90       include/inc_tracer.h
	cp include/inc_index_$(index).f90         include/inc_index.h
	cp include/inc_precision_$(precision).f90 include/inc_precision.h

makebin: $(PRJS)
	install $(PRJ1) $(BINDIR)/$(BINNAME)
	install $(PRJ2) $(BINDIR)/$(INITNAME)

$(PRJ1): $(PRJ1).o $(LIBDIR)/$(LIBNAME)
	$(LD) $(LFLAGS) -o $(PRJ1) $^

$(PRJ2): $(PRJ2).o $(LIBDIR)/$(LIBNAME)
	$(LD) $(LFLAGS) -o $(PRJ2) $^

$(PRJ1).o: $(PRJ1).f90 $(LIBNAME)
$(PRJ2).o: $(PRJ2).f90 $(LIBNAME)

$(LIBDIR)/$(LIBNAME): $(LIBNAME)
	mkdir -p $(LIBDIR)
	install $< $@
	$(RANLIB) $@

$(LIBNAME): $(MODS)
	$(AR) $@ $?

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
mod_monitor.o	: mod_monitor.f90 mod_stdio.o mod_fileio_h.o mod_process.o mod_geometrics.o mod_time.o
mod_atmos_sub_thermodyn.o	: mod_atmos_sub_thermodyn.f90 mod_stdio.o mod_const.o mod_time.o
mod_atmos_vars.o	: mod_atmos_vars.f90 mod_stdio.o mod_fileio_h.o mod_process.o mod_const.o mod_history.o mod_monitor.o mod_comm.o mod_fileio.o mod_time.o mod_grid_cartesian_plane.o mod_atmos_sub_thermodyn.o mod_atmos_sub_saturation.o
mod_atmos_sub_refstate.o	: mod_atmos_sub_refstate.f90 mod_stdio.o mod_process.o mod_fileio.o mod_fileio_h.o mod_const.o mod_grid_cartesian_plane.o
mod_atmos_sub_boundary.o	: mod_atmos_sub_boundary.f90 mod_stdio.o mod_fileio_h.o mod_process.o mod_const.o mod_grid_cartesian_plane.o mod_fileio.o mod_comm.o mod_time.o mod_atmos_sub_refstate.o
mod_atmos_dyn_fent_fullfct_lss.o	: mod_atmos_dyn_fent_fullfct_lss.f90 mod_stdio.o mod_time.o mod_process.o mod_const.o mod_grid_cartesian_plane.o mod_geometrics.o mod_comm.o mod_atmos_vars.o mod_atmos_sub_refstate.o mod_atmos_sub_boundary.o
mod_ocean_vars.o	: mod_ocean_vars.f90 mod_stdio.o mod_fileio_h.o mod_process.o mod_comm.o mod_fileio.o mod_time.o
mod_atmos_vars_sf.o	: mod_atmos_vars_sf.f90 mod_stdio.o mod_fileio_h.o mod_process.o mod_comm.o mod_fileio.o mod_time.o
mod_atmos_phy_sf_$(sfcflx).o	: mod_atmos_phy_sf_$(sfcflx).f90 mod_stdio.o mod_process.o mod_grid_cartesian_plane.o mod_const.o mod_atmos_vars.o mod_ocean_vars.o mod_atmos_vars_sf.o
mod_atmos_phy_tb_smg.o	: mod_atmos_phy_tb_smg.f90 mod_stdio.o mod_const.o mod_time.o mod_comm.o mod_grid_cartesian_plane.o mod_history.o mod_atmos_vars.o mod_atmos_vars_sf.o
mod_atmos_phy_mp_$(tracer).o	: mod_atmos_phy_mp_$(tracer).f90 mod_stdio.o mod_const.o mod_time.o mod_process.o mod_grid_cartesian_plane.o mod_history.o mod_atmos_vars.o mod_atmos_sub_precipitation.o mod_atmos_sub_thermodyn.o mod_atmos_sub_saturation.o mod_atmos_mp_saturation.o
mod_atmos_mp_saturation.o	: mod_atmos_mp_saturation.f90 mod_stdio.o mod_const.o mod_time.o
mod_atmos_phy_rd_$(radiation).o	: mod_atmos_phy_rd_$(radiation).f90 mod_stdio.o mod_process.o mod_comm.o mod_const.o mod_time.o mod_grid_cartesian_plane.o mod_history.o mod_atmos_vars.o 
mod_atmos.o	: mod_atmos.f90 mod_stdio.o mod_time.o mod_atmos_vars.o mod_atmos_sub_refstate.o mod_atmos_sub_boundary.o mod_atmos_sub_thermodyn.o mod_atmos_dyn_fent_fullfct_lss.o mod_atmos_phy_sf_$(sfcflx).o mod_atmos_phy_tb_smg.o mod_atmos_phy_mp_$(tracer).o mod_atmos_phy_rd_$(radiation).o
mod_ocean_fixedsst.o	: mod_ocean_fixedsst.f90 mod_stdio.o mod_process.o mod_ocean_vars.o mod_time.o mod_history.o
mod_ocean.o	: mod_ocean.f90 mod_stdio.o mod_time.o mod_ocean_vars.o mod_ocean_fixedsst.o
mod_atmos_sub_hydrostatic.o	: mod_atmos_sub_hydrostatic.f90 mod_stdio.o mod_const.o mod_comm.o mod_grid_cartesian_plane.o mod_process.o
mod_atmos_sub_saturation.o	: mod_atmos_sub_saturation.f90 mod_stdio.o mod_const.o mod_time.o
mod_atmos_sub_precipitation.o	: mod_atmos_sub_precipitation.f90 mod_stdio.o mod_time.o mod_const.o mod_grid_cartesian_plane.o mod_atmos_vars.o mod_atmos_sub_thermodyn.o
mod_mkinit.o	: mod_mkinit.f90 mod_stdio.o mod_process.o mod_const.o mod_random.o mod_grid_cartesian_plane.o mod_atmos_vars.o mod_atmos_sub_hydrostatic.o mod_atmos_sub_saturation.o

clean:
	rm -f $(PRJS) include/inc_index.h include/inc_tracer.h include/inc_precision.h *.mod *.o *.lst
distclean: clean
	rm -f $(LIBNAME)
	rm -f libscale_*.a

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.f90.o:
	$(FC) $(FFLAGS) -I./include -o $@ -c $<
.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

%.o: %.mod
