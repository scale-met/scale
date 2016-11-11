#! /bin/bash -x

GLEV=${1}
RLEV=${2}
TPROC=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=1"
#PJM --rsc-list "elapse=00:05:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin  "rank=* ${TOPDIR}/bin/${BINNAME} %r:./"
#PJM --stgin  "rank=* ./mkmnginfo.cnf          %r:./"
#PJM --stgout "rank=* %r:./*           ./"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
export XOS_MMM_L_ARENA_FREE=2

# run
mpiexec ./${BINNAME} mkmnginfo.cnf || exit

################################################################################
EOF1
