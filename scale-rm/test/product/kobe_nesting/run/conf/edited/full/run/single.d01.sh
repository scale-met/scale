#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=36"
#PJM --mpi "proc=36"
#PJM --rsc-list "elapse=00:30:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#-------------------------------------------------------------------------------
#PJM --stgin  "rank=* /data/ra000006/a00000/scale/bin/scale-rm        %r:./"
#PJM --stgin  "rank=* ../config/run.d01.s.conf                             %r:./"
#PJM --stgin  "rank=* ../init/domain_01/init_d01_00023587200.000.pe%06r.nc %r:./"
#PJM --stgin  "rank=* ../init/domain_01/boundary_d01.pe%06r.nc             %r:./"
#PJM --stgin  "rank=* ../input/domain_01/topo_d01.pe%06r.nc                %r:./"
#PJM --stgin  "rank=* ../input/domain_01/landuse_d01.pe%06r.nc             %r:./"
#-------------------------------------------------------------------------------
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/PARAG.29           %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/PARAPC.29          %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/VARDATA.RM29       %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/cira.nc            %r:./'
#PJM --stgin  'rank=* /data/ra000006/SCALE/database/rad/MIPAS/*      %r:./MIPAS/'
#-------------------------------------------------------------------------------
#PJM --stgout "rank=* %r:./*       ./domain_01/"
#PJM --stgout "rank=* %r:./prof/*  ./domain_01/prof/"
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
${fprof} mpiexec -n 36 ./scale-rm run.d01.s.conf || exit 1
