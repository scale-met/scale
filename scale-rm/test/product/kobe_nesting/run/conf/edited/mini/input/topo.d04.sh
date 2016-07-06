#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
#PJM --rsc-list "rscgrp=large"
#PJM --rsc-list "node=400"
#PJM --rsc-list "elapse=00:10:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#-------------------------------------------------------------------------------
#PJM --stgin  'rank=* /data/ra000006/a00000/scale/src/bin/scale-rm_pp        %r:./'
#PJM --stgin  'rank=* ../config/pp.d04.topo.conf                              %r:./'
#PJM --stgin  "rank=* /data/ra000006/SCALE/database/topo/DEM50M/Products/*    %r:./topo/"
#xxx --stgin  "rank=* /data/ra000006/SCALE/database/topo/GTOPO30/Products/*   %r:./topo/"
#-------------------------------------------------------------------------------
#PJM --stgout "rank=* %r:./* ./domain_04/"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
#
# run
mpiexec ./scale-rm_pp pp.d04.topo.conf   || exit 1
