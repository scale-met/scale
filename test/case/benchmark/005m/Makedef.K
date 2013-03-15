#
# ------ FOR K computer -----
#

##### for computation

FFLAGS_FAST  = -Kfast,parallel,ocl,preex,array_private,noalias=s,mfunc=2,dynamic_iteration \
               -Qi -Qt -X03 -Ncompdisp -Koptmsg=1 -Cpp \
               -Kprefetch_cache_level=all,prefetch_iteration_L2=50 -Ksimd

FFLAGS_DEBUG = -O0 \
               -Qi -Qt -X03 -v03s -v03d -v03o -Ncompdisp -Koptmsg=1 -Cpp \
               -Ec -Eg -Ha -He -Hf -Ho -Hs -Hu -Hx -Ncheck_global

#               -DDEBUG

# Performance monitor
# disable
#PERF_MONIT = -Ntl_notrt -U_FIPP_ -U_FAPP_
# fipp
PERF_MONIT = -Ntl_trt -D_FIPP_ -U_FAPP_
# fapp
#PERF_MONIT = -Ntl_trt -U_FIPP_ -D_FAPP_

# RDMA
# disable
RDMA = -U_USE_RDMA
# enable
#RDMA = -D_USE_RDMA

FC     = mpifrtpx
FFLAGS = $(FFLAGS_FAST) $(PERF_MONIT) $(RDMA)
#FFLAGS = $(FFLAGS_DEBUG) $(PERF_MONIT) $(RDMA)

FFLAGS_DYN ?= $(FFLAGS) -Knoprefetch,loop_nofission,loop_nofusion,nounroll -Ksimd=2

MODDIROPT ?= -M

CC     = mpifccpx
CFLAGS = -Kfast,parallel,ocl,preex,array_private,region_extension,restp=all -Qt -Ksimd $(PERF_MONIT)

LD     = $(FC)
LFLAGS = $(FFLAGS)

NETCDF_INCLUDE ?= -I/data/share005/cross-K/include
NETCDF_LIBS    ?= -L/data/share005/cross-K/lib -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz

##### for frontend
AR      = ar
ARFLAGS = r
RANLIB  = ranlib
JOBSUB  = pjsub
