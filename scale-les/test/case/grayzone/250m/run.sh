#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=16"
###PJM --rsc-list "node=1"
###PJM --rsc-list "elapse=00:05:00"
#PJM --rsc-list "elapse=24:00:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin  "rank=* ./scale-les_init %r:./"
#PJM --stgin  "rank=* ./scale-les  %r:./"
#PJM --stgin  "rank=* ../input_sounding.txt %r:./"
#PJM --stgin  "rank=* ../large_scale_w_force.txt %r:./"
#PJM --stgin  "rank=* ../sst_force.txt %r:./"
#PJM --stgin  "rank=*         ./init.conf %r:./"
#PJM --stgin  "rank=*         ./run.conf  %r:./"
#PJM --stgout "rank=* %r:./*      ./"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8

export LD_LIBRARY_PATH="/opt/aics/netcdf/k-serial/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY="/opt/aics/netcdf/k-serial/lib:$LD_LIBRARY"

# run
mpiexec ./scale-les_init init.conf || exit
ls -l 
echo "simulation starts"
mpiexec ./scale-les  run.conf  || exit
ls -l 

################################################################################
