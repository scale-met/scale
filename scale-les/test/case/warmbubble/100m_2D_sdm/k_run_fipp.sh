#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
##PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=1"
#PJM --mpi "proc=8"
#PJM --rsc-list "elapse=00:20:00"
#PJM --mail-list shinichiro.shima@gmail.com
#PJM --mpi "use-rankdir"
#PJM --stg-transfiles all
#PJM --stgin  "rank=* scale-les  %r:./"
#PJM --stgin  "rank=* run.conf  %r:./"
#PJM --stgin  "rank=* init_00000000000.000.pe%06r.nc %r:./"
#PJM --stgin  "rank=* random_number_init.pe%06r %r:./"
#PJM --stgout "rank=* %r:./* ./"
#PJM --stgout "rank=* %r:./Fprofd/* ./Fprofd/"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=1
export OMP_NUM_THREADS=1
export FLIB_CNTL_BARRIER_ERR="FALSE"

# run
# all range
fipp -C -d Fprofd -Icall,hwm mpiexec ./scale-les run.conf
# selected range
#fipp -C -Srange -d Fprofd -Icall,hwm mpiexec ./scale-les run.conf
#mpiexec ./scale-les run.conf

################################################################################
