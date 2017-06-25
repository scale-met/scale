################################################################################
#
# ------ For K computer -----
#
################################################################################

##### Fortran setting
#FFLAGS is set in Mkinclude. FFLAGS_DEBUG is used if SCALE_DEBUG=T

FC  = mpifrtpx
SFC = frtpx
MODDIROPT = -M

FFLAGS_FAST  = -Kfast,parallel,ocl,preex,array_private,noalias=s,mfunc=2 \
               -Qi -Qt -X03 -Ncompdisp -Koptmsg=1 -Cpp                   \
               -x-                                                       \
               -Ksimd                                                    \
               -Kauto,threadsafe

# NOTE: -Kauto, -Kthreadsafe, -Kautoobjstack, -Ktemparraystack depends on the stacksize.
#       If you find the segmentation fault error, please check.
FFLAGS_AGGRESSIVE = $(FFLAGS_FAST)  \
                   -Kautoobjstack   \
                   -Ktemparraystack \
                   -Kinstance=8     \
                   -Kswp_strong

FFLAGS_QUICKDEBUG = $(FFLAGS_FAST) -Nquickdbg -NRtrap

FFLAGS_DYN   = $(FFLAGS) -Knoprefetch,loop_nofission,loop_nofusion,nounroll,parallel_strong -Ksimd=2 -x100

FFLAGS_DEBUG = -O0                                                       \
               -Qi -Qt -X03 -v03s -v03d -v03o -Ncompdisp -Koptmsg=1 -Cpp \
               -Ec -Eg -Ha -He -Hf -Ho -Hs -Hu -Hx -Ncheck_global



##### C setting
#CFLAGS is set in Mkinclude. CFLAGS_DEBUG is used if SCALE_DEBUG=T

CC = mpifccpx

CFLAGS_FAST  = -Kfast,parallel,ocl,preex,array_private,region_extension,restp=all -Ksimd
CFLAGS_DEBUG = -O0



##### Special setting

ifeq ($(SCALE_ENABLE_OPENMP),T)
	OMP = -Kopenmp
# enable + disable parallel
#  OMP = -Kopenmp,noparallel
else
	OMP =
endif

### Performance monitor
# disable
PERF_MONIT = -Ntl_notrt -U_FIPP_ -U_FAPP_ -U_FINEPA_
# fipp
#PERF_MONIT = -Ntl_trt -D_FIPP_ -U_FAPP_ -U_FINEPA_
# fapp
#PERF_MONIT = -Ntl_trt -U_FIPP_ -D_FAPP_ -U_FINEPA_
# fine PA
#PERF_MONIT = -Ntl_trt -U_FIPP_ -U_FAPP_ -D_FINEPA_

# RDMA
# disable
RDMA = -U_USE_RDMA
# enable
#RDMA = -D_USE_RDMA

FFLAGS_SYSDEPEND = $(OMP) $(PERF_MONIT) $(RDMA)
CFLAGS_SYSDEPEND = $(OMP) $(PERF_MONIT) $(RDMA)



##### Linker setting

LD      = $(FC)
LDFLAGS = $(FFLAGS)



################################################################################

###### NetCDF library
# Default settings
SCALE_NETCDF_INCLUDE ?= -I/opt/aics/scale/include -I/opt/aics/netcdf/k-serial-noszip/include
SCALE_NETCDF_LIBS    ?= -L/opt/aics/scale/lib -L/opt/aics/netcdf/k-serial-noszip/lib-static -lpnetcdf -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm
SCALE_ENABLE_PNETCDF ?= T
# Or, please set environment variable
# example:
# export SCALE_NETCDF_INCLUDE="-I/opt/aics/netcdf/k-serial-noszip/include"
# export SCALE_NETCDF_LIBS="-L/opt/aics/netcdf/k-serial-noszip/lib-static -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm"

###### Math library
# Default settings
SCALE_MATHLIB_LIBS ?=
# Or, please set environment variable
# example:
# export SCALE_MATHLIB_LIBS="-SSL2"

##### for frontend
INSTALL = install
AR      = ar
ARFLAGS = r
RANLIB  = ranlib
JOBSUB  = pjsub
