#### Makefile ------------------------------------------------------------------------
#### 2014/11/26: Original (Ryuji Yoshida)


#### Intel Compiler with Intel MPI ---------------------------------------------------
FC=mpiifort

FFLAGS_CORE = -assume byterecl -convert big_endian -mcmodel=medium -heap-arrays

FFLAGS_FAST  = -m64 $(FFLAGS_CORE) -O3 -xHost -fp-model source

FFLAGS_DEBUG = -m64 $(FFLAGS_CORE) -O0 -fpe0 -fp-model source \
               -check all -warn all -g -traceback -assume nobuffered_io 
####----------------------------------------------------------------------------------

#### GNU Compiler with openMPI -------------------------------------------------------
#FC=mpif90
#FFLAGS_CORE = -m64 -std=f2003 -fall-intrinsics -pedantic-errors -fimplicit-none \
#              -fmodule-private -fconvert=big-endian -frecord-marker=4

#FFLAGS_FAST  = -O3 $(FFLAGS_CORE) -mtune=native -mcmodel=medium \
#               -funroll-loops -finline-functions -fprefetch-loop-arrays

#FFLAGS_DEBUG = -O0 $(FFLAGS_CORE) -ffpe-trap=invalid,zero,overflow  -fbounds-check \
#               -Wall -Wextra -Wcharacter-truncation -Wunderflow -g -fbacktrace
####----------------------------------------------------------------------------------

#### K computer & FX10 ---------------------------------------------------------------
#FC=mpifrtpx

#FFLAGS_CORE = -Qi -Qt -X03 -Ncompdisp -Koptmsg=1

#FFLAGS_FAST  = -Kfast,parallel,ocl,preex,array_private,noalias=s,mfunc=2 \
#               $(FFLAGS_CORE) -x- -Ksimd -Kauto,threadsafe

#FFLAGS_DEBUG = -O0 -v03s -v03d -v03o $(FFLAGS_CORE) \
#               -Ec -Eg -Ha -He -Hf -Ho -Hs -Hu -Hx -Ncheck_global
####----------------------------------------------------------------------------------

#FFLAGS= $(FFLAGS_FAST)
FFLAGS= $(FFLAGS_DEBUG)

LIBS=-L./
INCLUDES=-I./

OBJ=launcher
MODS=mod_vars.o mod_prof.o mod_comm.o mod_driver.o

####----------------------------------------------------------------------------------

all: launcher

launcher: launcher.o
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBS) -o $(OBJ) $(MODS) launcher.o

####----------------------------------------------------------------------------------

mod_vars.o: mod_vars.f90
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_prof.o: mod_prof.f90 mod_vars.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_comm.o: mod_comm.f90 mod_vars.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
mod_driver.o: mod_driver.f90 mod_vars.o mod_prof.o mod_comm.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
launcher.o: launcher.f90 mod_vars.o mod_prof.o mod_comm.o mod_driver.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ -c $<

####----------------------------------------------------------------------------------

clean:
	rm -f *.o *.mod $(OBJ)

####----------------------------------------------------------------------------------
