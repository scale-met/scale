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

if [ ${xy} -gt 1024 ]; then
   echo "Node usage is less than 1024. STOP"
   exit
fi

# Generate run.sh

cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
#PJM --rsc-list "rscgrp=micro"
#PJM --rsc-list "node=${TPROC}"
#PJM --rsc-list "elapse=00:25:00"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
#export fu08bf=1

rm -f             ./input
ln -svf ${DATDIR} ./input

# run
EOF1

if [ ! ${PPNAME}   = "NONE" ]; then
   echo "${MPIEXEC} ${BINDIR}/${PPNAME}   ${PPCONF}   || exit 1" >> ./run.sh
fi

if [ ! ${INITNAME} = "NONE" ]; then
   echo "${MPIEXEC} ${BINDIR}/${INITNAME} ${INITCONF} || exit 1" >> ./run.sh
fi

if [ ! ${BINNAME} = "NONE" ]; then
   echo "${MPIEXEC} ${BINDIR}/${BINNAME}  ${RUNCONF}  || exit 1" >> ./run.sh
fi
