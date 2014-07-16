#! /bin/bash -x
################################################################################
#
# ------ FOR Linux64 & intel C&fortran & intel mpi -----
#
################################################################################
export FORT_FMT_RECL=400


ln -svf ../../../data/rad/PARAG.29 .
ln -svf ../../../data/rad/PARAPC.29 .
ln -svf ../../../data/rad/VARDATA.RM29 .
ln -svf ../../../data/rad/cira.nc .
ln -svf ../../../data/rad/MIPAS/day.atm .
ln -svf ../../../data/rad/MIPAS/equ.atm .
ln -svf ../../../data/rad/MIPAS/sum.atm .
ln -svf ../../../data/rad/MIPAS/win.atm .

# run
mpirun -np 1 ./scale-les_init init.conf || exit
mpirun -np 1 ./scale-les  run.restart0.conf  || exit
mpirun -np 1 ./scale-les  run.restartA.conf  || exit
mpirun -np 1 ./scale-les  run.restartB.conf  || exit

# check
ruby diff.rb restart0_*.nc restartB_*.nc

################################################################################
