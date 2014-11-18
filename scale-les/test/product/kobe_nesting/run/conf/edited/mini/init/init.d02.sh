#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=4"
#PJM --rsc-list "elapse=00:30:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#-------------------------------------------------------------------------------
#PJM --stgin  "rank=* /data/ra000006/a00000/scale/src/bin/scale-les_init  %r:./"
#PJM --stgin  "rank=* ../config/init.d02.conf                             %r:./"
#PJM --stgin  "rank=* ./extdata/*.nc                                      %r:./"
#PJM --stgin  "rank=* ../input/domain_02/topo_d02.pe%06r.nc               %r:./"
#PJM --stgin  "rank=* ../input/domain_02/landuse_d02.pe%06r.nc            %r:./"
#-------------------------------------------------------------------------------
#PJM --stgout "rank=* %r:./*  ./domain_02/"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8

# run
mpiexec ./scale-les_init init.d02.conf || exit 1
