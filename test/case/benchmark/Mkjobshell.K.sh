#! /bin/bash -x

# Arguments
BINDIR=${1}
INITNAME=${2}
BINNAME=${3}
INITCONF=${4}
RUNCONF=${5}
TPROC=${6}

# System specific
MPIEXEC="mpiexec"

array=( `echo ${TPROC} | tr -s 'x' ' '`)
x=${array[0]}
y=${array[1]:-1}
let xy="${x} * ${y}"

if [ ${xy} -gt 36864 ]; then
   rscgrp="huge"
elif [ ${xy} -gt 384 ]; then
   rscgrp="large"
else
   rscgrp="small"
fi

# Generate run.sh

cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
#PJM --rsc-list "rscgrp=${rscgrp}"
#PJM --rsc-list "node=${TPROC}"
#PJM --rsc-list "elapse=01:00:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin  "rank=* ${BINDIR}/${INITNAME} %r:./"
#PJM --stgin  "rank=* ${BINDIR}/${BINNAME}  %r:./"
#PJM --stgin  "rank=*         ./${INITCONF} %r:./"
#PJM --stgin  "rank=*         ./${RUNCONF}  %r:./"
#PJM --stgout "rank=* %r:./*      ./"
#PJM --stgout "rank=* %r:./prof/* ./prof/"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
#export fu08bf=1

rm -rf ./prof

# run
                              ${MPIEXEC} ./${INITNAME} ${INITCONF} || exit
fipp -C -Srange -Ihwm -d prof ${MPIEXEC} ./${BINNAME}  ${RUNCONF}  || exit

################################################################################
EOF1
