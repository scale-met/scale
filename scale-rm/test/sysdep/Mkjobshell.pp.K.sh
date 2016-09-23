#! /bin/bash -x

# Arguments
BINDIR=${1}
PPNAME=${2}
INITNAME=${3}
BINNAME=${4}
PPCONF=${5}
INITCONF=${6}
RUNCONF=${7}
TPROC=${8}
DATDIR=${9}
DATPARAM=(`echo ${10} | tr -s ',' ' '`)
DATDISTS=(`echo ${11} | tr -s ',' ' '`)

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
#PJM --rsc-list "elapse=02:00:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
EOF1

if [ ! ${PPNAME}   = "NONE" ]; then
   echo "#PJM --stgin  'rank=* ${BINDIR}/${PPNAME}   %r:./'" >> ./run.sh
   echo "#PJM --stgin  'rank=*         ./${PPCONF}   %r:./'" >> ./run.sh
fi
if [ ! ${INITNAME}   = "NONE" ]; then
   echo "#PJM --stgin  'rank=* ${BINDIR}/${INITNAME} %r:./'" >> ./run.sh
   echo "#PJM --stgin  'rank=*         ./${INITCONF} %r:./'" >> ./run.sh
fi
if [ ! ${BINNAME}   = "NONE" ]; then
   echo "#PJM --stgin  'rank=* ${BINDIR}/${BINNAME}  %r:./'" >> ./run.sh
   echo "#PJM --stgin  'rank=*         ./${RUNCONF}  %r:./'" >> ./run.sh
fi

cat << EOF2 >> ./run.sh
#PJM --stgin  "rank=* ${DATDIR}/* %r:./input/"
#PJM --stgout "rank=* %r:./* ./"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
#export fu08bf=1

# run
EOF2

if [ ! ${PPNAME}   = "NONE" ]; then
   echo "${MPIEXEC} ./${PPNAME}   ${PPCONF}   || exit 1" >> ./run.sh
fi

if [ ! ${INITNAME} = "NONE" ]; then
   echo "${MPIEXEC} ./${INITNAME} ${INITCONF} || exit 1" >> ./run.sh
fi

if [ ! ${BINNAME} = "NONE" ]; then
   echo "${MPIEXEC} ./${BINNAME}  ${RUNCONF}  || exit 1" >> ./run.sh
fi
