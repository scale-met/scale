#! /bin/bash -x
################################################################################
#
# ------ FOR Linux64 & gnu C&fortran & openmpi -----
#
################################################################################
export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y


# run
mpirun -np  2 /home/tyamaura/scale/scale-les/test/case/coldbubble/sp_1000m_250m/scale-les_init init.d01.conf || exit
mpirun -np  4 /home/tyamaura/scale/scale-les/test/case/coldbubble/sp_1000m_250m/scale-les_init init.d02.conf || exit
mpirun -np  6 /home/tyamaura/scale/scale-les/test/case/coldbubble/sp_1000m_250m/scale-les launch.conf || exit

################################################################################
