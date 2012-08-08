################################################################################
#
# Makefile for main program
#
################################################################################

include ./sysdep/Makedef.$(SCALE_SYS)

VPATH = include:common:io:grid:atmos:ocean:preprocess

LIBDIR = ../lib

#default index & tracer set
index   = 1256x8x8
tracer  = kessler
precision = real8
postfix = 

PRJ1 = scale3
PRJ2 = scale3_init

PRJS = $(PRJ1) $(PRJ2)

LIBNAME = libscale_$(index)_$(tracer)_$(precision)_$(postfix).a

MODS =	\
	mod_stdio.o	\
	mod_process.o	\
	mod_const.o	\
	mod_random.o	\
	mod_time.o	\
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
	mod_atmos_dyn_fent_fct.o	\
	mod_atmos_phy_sf.o	\
	mod_atmos_vars_sf.o	\
	mod_atmos_phy_tb_smg.o	\
	mod_atmos_phy_mp_kes.o  \
	mod_ocean.o	\
	mod_ocean_vars.o	\
	mod_ocean_fixedsst.o	\
	mod_mkinit.o	\
	mod_atmos_sub_hydrostatic.o	\
	mod_atmos_sub_saturation.o

all: makegtool makedir makeinclude makebin

makegtool:
	cd gtool; make

makedir:
	mkdir -p ../bin/$(SCALE_SYS)

makeinclude:
	cp include/inc_tracer_$(tracer).f90 include/inc_tracer.h
	cp include/inc_index_$(index).f90   include/inc_index.h
	cp include/inc_precision_$(precision).f90 include/inc_precision.h

makebin: $(PRJS)
	install $(PRJ1) ../bin/$(SCALE_SYS)/$(PRJ1)_$(index)_$(tracer)_$(precision)_$(postfix)
	install $(PRJ2) ../bin/$(SCALE_SYS)/$(PRJ2)_$(index)_$(tracer)_$(precision)_$(postfix)

$(PRJ1): $(PRJ1).o $(LIBDIR)/$(LIBNAME) $(LIBDIR)/libgtool.a
	$(LD) $(LFLAGS) -o $(PRJ1) $^ $(NETCDF_LIBS)

$(PRJ2): $(PRJ2).o $(LIBDIR)/$(LIBNAME) $(LIBDIR)/libgtool.a
	$(LD) $(LFLAGS) -o $(PRJ2) $^ $(NETCDF_LIBS)

$(PRJ1).o: $(PRJ1).f90 libscale.a
$(PRJ2).o: $(PRJ2).f90 libscale.a

$(LIBDIR)/libgtool.a: makegtool

$(LIBDIR)/$(LIBNAME): libscale.a
	mkdir -p $(LIBDIR)
	install libscale.a $(LIBDIR)/$(LIBNAME)

libscale.a: $(MODS)
	ar r libscale.a $?

mod_stdio.o	: mod_stdio.f90
mod_process.o	: mod_process.f90 mod_stdio.o
mod_const.o	: mod_const.f90 mod_stdio.o mod_process.o
mod_random.o	: mod_random.f90 mod_stdio.o mod_process.o
mod_time.o	: mod_time.f90 mod_stdio.o mod_process.o
mod_grid_cartesian_plane.o	: mod_grid_cartesian_plane.f90 mod_stdio.o mod_process.o $(LIBDIR)/libgtool.a
mod_geometrics.o	: mod_geometrics.f90 mod_stdio.o mod_process.o mod_time.o mod_grid_cartesian_plane.o mod_const.o
mod_comm.o	: mod_comm.f90 mod_stdio.o mod_time.o mod_process.o mod_const.o mod_geometrics.o
mod_history.o	: mod_history.f90 mod_stdio.o $(LIBDIR)/libgtool.a
mod_monitor.o	: mod_monitor.f90 mod_stdio.o mod_process.o mod_geometrics.o mod_time.o $(LIBDIR)/libgtool.a
mod_atmos_sub_thermodyn.o	: mod_atmos_sub_thermodyn.f90 mod_stdio.o mod_const.o mod_time.o
mod_atmos_vars.o	: mod_atmos_vars.f90 mod_stdio.o mod_process.o mod_const.o mod_history.o mod_monitor.o mod_comm.o mod_time.o mod_grid_cartesian_plane.o mod_atmos_sub_thermodyn.o mod_atmos_sub_saturation.o
mod_atmos_sub_refstate.o	: mod_atmos_sub_refstate.f90 mod_stdio.o mod_process.o mod_const.o mod_grid_cartesian_plane.o
mod_atmos_sub_boundary.o	: mod_atmos_sub_boundary.f90 mod_stdio.o mod_process.o mod_const.o mod_grid_cartesian_plane.o mod_comm.o mod_time.o mod_atmos_sub_refstate.o
mod_atmos_dyn_fent_fct.o	: mod_atmos_dyn_fent_fct.f90 mod_stdio.o mod_time.o mod_process.o mod_const.o mod_grid_cartesian_plane.o mod_geometrics.o mod_comm.o mod_atmos_vars.o mod_atmos_sub_refstate.o mod_atmos_sub_boundary.o
mod_ocean_vars.o	: mod_ocean_vars.f90 mod_stdio.o mod_process.o mod_comm.o mod_time.o
mod_atmos_vars_sf.o	: mod_atmos_vars_sf.f90 mod_stdio.o mod_process.o mod_comm.o mod_time.o
mod_atmos_phy_sf.o	: mod_atmos_phy_sf.f90 mod_stdio.o mod_process.o mod_grid_cartesian_plane.o mod_const.o mod_atmos_vars.o mod_ocean_vars.o mod_atmos_vars_sf.o
mod_atmos_phy_tb_smg.o	: mod_atmos_phy_tb_smg.f90 mod_stdio.o mod_const.o mod_time.o mod_comm.o mod_grid_cartesian_plane.o mod_history.o mod_atmos_vars.o mod_atmos_vars_sf.o
mod_atmos_phy_mp_kes.o  : mod_atmos_phy_mp_kes.f90 mod_stdio.o mod_process.o mod_const.o mod_time.o mod_grid_cartesian_plane.o mod_history.o mod_atmos_vars.o
mod_atmos.o	: mod_atmos.f90 mod_stdio.o mod_time.o mod_atmos_vars.o mod_atmos_sub_refstate.o mod_atmos_sub_boundary.o mod_atmos_sub_thermodyn.o mod_atmos_dyn_fent_fct.o mod_atmos_phy_sf.o mod_atmos_phy_tb_smg.o mod_atmos_phy_mp_kes.o
mod_ocean_fixedsst.o	: mod_ocean_fixedsst.f90 mod_stdio.o mod_process.o mod_ocean_vars.o mod_time.o mod_history.o
mod_ocean.o	: mod_ocean.f90 mod_stdio.o mod_time.o mod_ocean_vars.o mod_ocean_fixedsst.o
mod_atmos_sub_hydrostatic.o	: mod_atmos_sub_hydrostatic.f90 mod_stdio.o mod_const.o mod_comm.o mod_grid_cartesian_plane.o mod_process.o
mod_atmos_sub_saturation.o	: mod_atmos_sub_saturation.f90 mod_stdio.o mod_const.o mod_time.o
mod_mkinit.o	: mod_mkinit.f90 mod_stdio.o mod_process.o mod_const.o mod_random.o mod_grid_cartesian_plane.o mod_atmos_vars.o mod_atmos_sub_hydrostatic.o mod_atmos_sub_saturation.o

clean:
	rm -f $(MODS) *.mod

distclean: clean
	rm -f $(PRJS) include/inc_index.h include/inc_tracer.h include/inc_precision.h libscale.a

.SUFFIXES:
.SUFFIXES: .o .f90 .c

.f90.o:
	$(FC) $(FFLAGS) -I./include -I./gtool -o $@ -c $<
.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

%.o: %.mod


.PHONY : clean distclean
