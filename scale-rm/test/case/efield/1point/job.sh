#! /bin/bash -x
################################################################################
#
# ------ For Linux64 & gnu fortran&C & openmpi -----
#
################################################################################
#export FORT_FMT_RECL=400
#export GFORTRAN_UNBUFFERED_ALL=Y

# run
mpirun -np 6 --oversubscribe ./scale-rm_init init.conf || exit 1
mpirun -np 6 --oversubscribe ./scale-rm run.conf || exit 1


################################################################################
