#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
#PJM --rsc-list "rscgrp=large"
#PJM --rsc-list "node=14400"
#PJM --mpi "proc=14400"
#PJM --rsc-list "elapse=00:30:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#-------------------------------------------------------------------------------
#PJM --stgin  "rank=* /data/ra000006/a00000/scale/src/bin/scale-les        %r:./"
#PJM --stgin  "rank=* ../config/run.d04.s.conf                             %r:./"
#PJM --stgin  "rank=* ../init/domain_04/init_d04_00023587200.000.pe%06r.nc %r:./"
#PJM --stgin  "rank=* ../init/domain_04/boundary_d04.pe%06r.nc             %r:./"
#PJM --stgin  "rank=* ../input/domain_04/topo_d04.pe%06r.nc                %r:./"
#PJM --stgin  "rank=* ../input/domain_04/landuse_d04.pe%06r.nc             %r:./"
#-------------------------------------------------------------------------------
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/PARAG.29           %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/PARAPC.29          %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/VARDATA.RM29       %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/cira.nc            %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/MIPAS/*      %r:./MIPAS/'
#-------------------------------------------------------------------------------
#PJM --stgout "rank=* %r:./*       ./domain_04/"
#PJM --stgout "rank=* %r:./prof/*  ./domain_04/prof/"
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
${fprof} mpiexec -n 14400 ./scale-les run.d04.s.conf || exit 1
