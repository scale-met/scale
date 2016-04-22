#! /bin/bash -x
################################################################################
#
# ------ FOR Linux64 & intel C&fortran & intel mpi -----
#
################################################################################
export FORT_FMT_RECL=400


#rm -f             ./input
#ln -svf /data0/scale_database/topo/GTOPO30/Products ./input

# run
mpirun -np 36 ./scale-les_pp   pp.d01.conf   || exit 1
mpirun -np 36 ./scale-les_pp   pp.d02.conf   || exit 1
mpirun -np 48 ./scale-les_pp   pp.d03.conf   || exit 1
