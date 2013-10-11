#! /bin/bash -x
################################################################################
#
# ------ FOR MacOSX & gfortran4.6 & OpenMPI1.6 -----
#
################################################################################
export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

ln -sv ../../../../bin/MacOSX-gnu-ompi/scale3_init_200x8x8_real8_heve_sn13_const_mstrnX_smg_ .
ln -sv ../../../../bin/MacOSX-gnu-ompi/scale3_200x8x8_real8_heve_sn13_const_mstrnX_smg_  .

# run
openmpirun -np 6 ./scale3_init_200x8x8_real8_heve_sn13_const_mstrnX_smg_ init.conf || exit
openmpirun -np 6 ./scale3_200x8x8_real8_heve_sn13_const_mstrnX_smg_  run.conf  || exit

################################################################################
