#! /bin/bash -x
################################################################################
#
# ------ FOR Linux64 & GNU C&fortran & openmpi -----
#
################################################################################
export FORT_FMT_RECL=400


rm -f             ./input
ln -svf /data3/kenshi/scale_database/topo/GTOPO30/Products ./input

# run
mpirun -np 4 /home/tyamaura/scale/scale-rm/test/case_real/check_mass/topo/scale-rm_pp   pp.d01.conf   || exit 1
mpirun -np 4 /home/tyamaura/scale/scale-rm/test/case_real/check_mass/topo/scale-rm_pp   pp.d02.conf   || exit 1
