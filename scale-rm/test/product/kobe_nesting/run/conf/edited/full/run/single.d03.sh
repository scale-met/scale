#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
#PJM --rsc-list "rscgrp=large"
#PJM --rsc-list "node=576"
#PJM --mpi "proc=576"
#PJM --rsc-list "elapse=00:30:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#-------------------------------------------------------------------------------
#PJM --stgin  "rank=* /data/ra000006/a00000/scale/src/bin/scale-rm        %r:./"
#PJM --stgin  "rank=* ../config/run.d03.s.conf                             %r:./"
#PJM --stgin  "rank=* ../init/domain_03/init_d03_00023587200.000.pe%06r.nc %r:./"
#PJM --stgin  "rank=* ../init/domain_03/boundary_d03.pe%06r.nc             %r:./"
#PJM --stgin  "rank=* ../input/domain_03/topo_d03.pe%06r.nc                %r:./"
#PJM --stgin  "rank=* ../input/domain_03/landuse_d03.pe%06r.nc             %r:./"
#-------------------------------------------------------------------------------
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/PARAG.29           %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/PARAPC.29          %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/VARDATA.RM29       %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/cira.nc            %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/MIPAS/*      %r:./MIPAS/'
#-------------------------------------------------------------------------------
#PJM --stgout "rank=* %r:./*       ./domain_03/"
#PJM --stgout "rank=* %r:./prof/*  ./domain_03/prof/"
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
${fprof} mpiexec -n 576 ./scale-rm run.d03.s.conf || exit 1
