# ------ FOR K -----

FC     = mpifrtpx
# new setting 20120413
FFLAGS = -Kfast,parallel,ocl,preex,array_private,auto,noalias=s,parallel_compress,parallel_iteration=8,instance=8 \
         -Qi -Qt -X03 -v03s -v03d -v03o -Ncompdisp -Koptmsg=1 -Cpp \
         -Kprefetch_cache_level=all,prefetch_iteration_L2=50 -Ksimd -Kdynamic_iteration

FFLAGS_DYN = -Kfast,parallel,ocl,preex,array_private,auto,noalias=s,parallel_compress,parallel_iteration=8,instance=8 \
             -Qi -Qt -X03 -v03s -v03d -v03o -Ncompdisp -Koptmsg=1 -Cpp \
             -Knoprefetch -Ksimd=2 -Kdynamic_iteration

# old setting 20120214
#FFLAGS = -Kvisimpact,auto,openmp,array_private,ilfunc,ocl,mfunc=2,preex,noalias=s,parallel_compress \
#         -Ntl_trt -Qt -V -Cpp -Ksimd -D_FPCOLL_ \
#         -Kprefetch_cache_level=all,prefetch_iteration_L2=50

#FFLAGS_DYN = -Kvisimpact,auto,openmp,array_private,ilfunc,ocl,mfunc=2,preex,noalias=s,parallel_compress \
#             -Knoprefetch,loop_nofission,loop_nofusion,nounroll \
#             -Ntl_trt -Qt -V -Cpp -Ksimd=2 -D_FPCOLL_

# old setting 20120104
#FFLAGS = -Kvisimpact,ocl,ilfunc,preex,array_private,prefetch_cache_level=all,auto,mfunc=2,noalias=s,instance=8,prefetch_iteration_L2=50 \
#         -Ntl_trt -Qt -V -Cpp -Ksimd -D_FPCOLL_ \
#         -v03s -Ncompdisp,maxserious=2,nofreealloc -Kautoobjstack

CC     = mpifccpx
CFLAGS = -Kfast,parallel,ocl,region_extension,array_private,preex,openmp,restp=all -Koptmsg=2 -Ksimd

LD     = $(FC)
LFLAGS = $(FFLAGS)

# netcdf
NETCDF_INCLUDE = -I/data/group084/netcdf/include
NETCDF_LIBS = -L/data/group084/netcdf/lib -lnetcdf -L/data/group084/hdf5.mpi/lib -lhdf5_hl -lhdf5 -lsz -lz
# for frontend
FRONT_CC     =  LD_LIBRARY_PATH=/home/group084/user0202/gnu/lib:/home/group084/user0202/gnu/lib64:$$LD_LIBRARY_PATH /home/group084/user0202/gnu/bin/gcc
FRONT_CFLAGS = -O3

FRONT_F90    = LD_LIBRARY_PATH=/home/group084/user0202/gnu/lib:/home/group084/user0202/gnu/lib64:$$LD_LIBRARY_PATH /home/group084/user0202/gnu/bin/gfortran
FRONT_FFLAGS = -O3

FRONT_LD     = LD_LIBRARY_PATH=/home/group084/user0202/gnu/lib:/home/group084/user0202/gnu/lib64:$$LD_LIBRARY_PATH /home/group084/user0202/gnu/bin/gfortran
FRONT_LFLAGS = -lm

AR     = ar r
RANLIB = ranlib
