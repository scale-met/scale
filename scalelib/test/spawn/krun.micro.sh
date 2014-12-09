#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
#PJM --rsc-list "rscgrp=micro"
#PJM --rsc-list "node=15"
#PJM --mpi "shape=2"
#PJM --mpi "proc=2"
#PJM --rsc-list "elapse=00:00:30"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
#export THREAD_STACK_SIZE=256000
#export fu08bf=1

fprof="fipp -C -Stotal -Icall,hwm -d prof"

# run

${fprof} mpiexec -n 2  ./driver run.d01.conf || exit 1

