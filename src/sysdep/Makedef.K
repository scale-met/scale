# ------ FOR K -----
#
FPCOLL=-U_FPCOLL_
SIMD=-Ksimd


FC     = mpifrtpx
FFLAGS = -Kvisimpact,prefetch_iteration_L2=50,prefetch_cache_level=all,auto,openmp,array_private,ilfunc,ocl,mfunc=2,preex,noalias=s,parallel_compress,parallel_strong \
         -Ntl_trt $(FPCOLL) $(SIMD) $(FJTIMER)
#
CC     = mpifccpx
CFLAGS = -Kvisimpact,openmp,array_private,ilfunc,ocl,mfunc,preex,restp=all -Ntl_trt 
#
#
LD = $(FC)
LFLAGS = -Kvisimpact,openmp,array_private,ilfunc,ocl,mfunc,preex -Ntl_trt

# for frontend
FRONT_CC     =  LD_LIBRARY_PATH=/home/user0114/gnu/lib:/home/user0114/gnu/lib64:$$LD_LIBRARY_PATH /home/user0114/gnu/bin/gcc
FRONT_CFLAGS = -O3

FRONT_F90    = LD_LIBRARY_PATH=/home/user0114/gnu/lib:/home/user0114/gnu/lib64:$$LD_LIBRARY_PATH /home/user0114/gnu/bin/gfortran
FRONT_FFLAGS = -O3

FRONT_LD     = LD_LIBRARY_PATH=/home/user0114/gnu/lib:/home/user0114/gnu/lib64:$$LD_LIBRARY_PATH /home/user0114/gnu/bin/gfortran
FRONT_LFLAGS = -lm
