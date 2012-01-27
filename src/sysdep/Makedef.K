# ------ FOR K -----

FC     = mpifrtpx
#FFLAGS = -Kvisimpact,ocl,ilfunc,preex,array_private,prefetch_cache_level=all,auto,mfunc=2,noalias=s,instance=8,prefetch_iteration_L2=50 \
#         -Ntl_trt -Qt -V -Cpp -Ksimd -D_FPCOLL_ \
#         -v03s -Ncompdisp,maxserious=2,nofreealloc -Kautoobjstack
FFLAGS = -Kvisimpact,prefetch_iteration_L2=50,prefetch_cache_level=all,auto,openmp,array_private,ilfunc,ocl,mfunc=2,preex,noalias=s,parallel_compress,parallel_strong \
         -Ntl_trt  -Qt -V -Ksimd -Cpp -D_FPCOLL_

CC     = mpifccpx
#CFLAGS = -Kfast,parallel,ocl,preex,array_private,auto -Qt -V -Ksimd=2
CFLAGS = -Kvisimpact,openmp,array_private,ilfunc,ocl,mfunc,preex,restp=all -Ntl_trt -Qt -V -Ksimd

LD     = $(FC)
LFLAGS = $(FFLAGS)

# for frontend
FRONT_CC     =  LD_LIBRARY_PATH=/home/user0114/gnu/lib:/home/user0114/gnu/lib64:$$LD_LIBRARY_PATH /home/user0114/gnu/bin/gcc
FRONT_CFLAGS = -O3

FRONT_F90    = LD_LIBRARY_PATH=/home/user0114/gnu/lib:/home/user0114/gnu/lib64:$$LD_LIBRARY_PATH /home/user0114/gnu/bin/gfortran
FRONT_FFLAGS = -O3

FRONT_LD     = LD_LIBRARY_PATH=/home/user0114/gnu/lib:/home/user0114/gnu/lib64:$$LD_LIBRARY_PATH /home/user0114/gnu/bin/gfortran
FRONT_LFLAGS = -lm
