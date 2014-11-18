#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=100"
#PJM --mpi "proc=100"
#PJM --rsc-list "elapse=00:30:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#-------------------------------------------------------------------------------
#PJM --stgin  "rank=* /data/ra000006/a00000/scale/src/bin/scale-les        %r:./"
#PJM --stgin  "rank=* ../config/run.d02.s.conf                             %r:./"
#PJM --stgin  "rank=* ../init/domain_02/init_d02_00023587200.000.pe%06r.nc %r:./"
#PJM --stgin  "rank=* ../init/domain_02/boundary_d02.pe%06r.nc             %r:./"
#PJM --stgin  "rank=* ../input/domain_02/topo_d02.pe%06r.nc                %r:./"
#PJM --stgin  "rank=* ../input/domain_02/landuse_d02.pe%06r.nc             %r:./"
#-------------------------------------------------------------------------------
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/PARAG.29           %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/PARAPC.29          %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/VARDATA.RM29       %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/cira.nc            %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/MIPAS/*      %r:./MIPAS/'
#-------------------------------------------------------------------------------
#PJM --stgout "rank=* %r:./*       ./domain_02/"
#PJM --stgout "rank=* %r:./prof/*  ./domain_02/prof/"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
#ls -lha
#ls -lha ./MIPAS/

fprof="fipp -C -Srange -Icall,hwm -d prof"

# run
${fprof} mpiexec -n 100 ./scale-les run.d02.s.conf || exit 1
