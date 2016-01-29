#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
##PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=1"
#PJM --mpi "proc=8"
#PJM --rsc-list "elapse=0:10:00"
#PJM --mpi "use-rankdir"
#PJM --stgin  "rank=* ./scale-les_init  %r:./"
#PJM --stgin  "rank=* ./init.conf  %r:./"
#PJM --stgout "rank=* %r:./* ./"
#PJM --mail-list shinichiro.shima@gmail.com
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=1
export OMP_NUM_THREADS=1
export FLIB_CNTL_BARRIER_ERR="FALSE"
#export LD_LIBRARY_PATH=/data/hp140094/k01949/applications/from_chihaya/netcdf-cross/Linux-s64fx/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/home/usr7/i70107a/applications/netcdf-cross/Linux-s64fx/lib:$LD_LIBRARY_PATH

# run
mpiexec ./scale-les_init init.conf

################################################################################
