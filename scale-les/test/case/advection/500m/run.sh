#! /bin/bash -x
################################################################################
#
# ------ FOR Linux64 & gnu C&fortran & openmpi -----
#
################################################################################
export FORT_FMT_RECL=400


ln -sv ../../../../bin/Linux64-gnu-ompi/scale3_init_500m40x40x40_real8_heve_sn13_const_dummy_dummy_ .
ln -sv ../../../../bin/Linux64-gnu-ompi/scale3_500m40x40x40_real8_heve_sn13_const_dummy_dummy_  .

# run
mpirun -np 6 ./scale3_init_500m40x40x40_real8_heve_sn13_const_dummy_dummy_ init.conf || exit
mpirun -np 6 ./scale3_500m40x40x40_real8_heve_sn13_const_dummy_dummy_  run.conf  || exit

################################################################################
