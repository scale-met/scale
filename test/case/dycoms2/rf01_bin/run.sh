#! /bin/bash -x
################################################################################
#
# ------ FOR Linux64 & intel C&fortran & intel mpi -----
#
################################################################################
export FORT_FMT_RECL=400


ln -sv ../../../../bin/Linux64-intel-impi/scale3_init_200x8x8_real8_heve_binw_const_dycoms2_smg_ .
ln -sv ../../../../bin/Linux64-intel-impi/scale3_200x8x8_real8_heve_binw_const_dycoms2_smg_  .

# run
mpirun -np 1 ./scale3_init_200x8x8_real8_heve_binw_const_dycoms2_smg_ init.conf || exit
mpirun -np 1 ./scale3_200x8x8_real8_heve_binw_const_dycoms2_smg_  run.conf  || exit

################################################################################
