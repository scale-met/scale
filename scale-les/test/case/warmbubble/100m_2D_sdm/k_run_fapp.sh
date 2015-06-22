#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
##PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=1"
#PJM --mpi "proc=8"
#PJM --rsc-list "elapse=01:00:00"
#PJM --mail-list shinichiro.shima@gmail.com
#PJM --mpi "use-rankdir"
#PJM --stg-transfiles all
#PJM --stgin  "rank=* scale-les  %r:./"
#PJM --stgin  "rank=* run.conf  %r:./"
#PJM --stgin  "rank=* init_00000000000.000.pe%06r.nc %r:./"
#PJM --stgin  "rank=* random_number_init.pe%06r %r:./"
#PJM --stgout "rank=* %r:./* ./"
#PJM --stgout "rank=* %r:./Fprofd_cache/* ./Fprofd_cache/"
#PJM --stgout "rank=* %r:./Fprofd_instr/* ./Fprofd_instr/"
#PJM --stgout "rank=* %r:./Fprofd_memac/* ./Fprofd_memac/"
#PJM --stgout "rank=* %r:./Fprofd_perfo/* ./Fprofd_perfo/"
#PJM --stgout "rank=* %r:./Fprofd_stati/* ./Fprofd_stati/"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=1
export OMP_NUM_THREADS=1
export FLIB_CNTL_BARRIER_ERR="FALSE"

# run
LD="./scale-les run.conf"
MPIEXEC="mpiexec"
#
fapp -C -d Fprofd_cache -L1 -Hevent=Cache ${MPIEXEC} ${LD}
fapp -C -d Fprofd_instr -L1 -Hevent=Instructions ${MPIEXEC} ${LD}
fapp -C -d Fprofd_memac -L1 -Hevent=MEM_access ${MPIEXEC} ${LD}
fapp -C -d Fprofd_perfo -L1 -Hevent=Performance ${MPIEXEC} ${LD}
fapp -C -d Fprofd_stati -L1 -Hevent=Statistics ${MPIEXEC} ${LD}
#mpiexec ./scale-les run.conf

################################################################################
