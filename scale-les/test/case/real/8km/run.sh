#! /bin/bash -x
################################################################################
#
# ------ FOR Linux64 & intel C&fortran & intel mpi -----
#
################################################################################
export FORT_FMT_RECL=400


# run
mpirun -np 16 /work2/scale/git/scale/scale-les/test/case/real/10km/scale-les_init init.conf || exit
mpirun -np 16 /work2/scale/git/scale/scale-les/test/case/real/10km/scale-les  run.conf  || exit

################################################################################
