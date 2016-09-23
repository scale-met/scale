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

# for Oakleaf-FX
# if [ ${xy} -gt 480 ]; then
#    rscgrp="x-large"
# elif [ ${xy} -gt 372 ]; then
#    rscgrp="large"
# elif [ ${xy} -gt 216 ]; then
#    rscgrp="medium"
# elif [ ${xy} -gt 12 ]; then
#    rscgrp="small"
# else
#    rscgrp="short"
# fi

# for AICS-FX10
if [ ${xy} -gt 96 ]; then
   rscgrp="huge"
elif [ ${xy} -gt 24 ]; then
   rscgrp="large"
else
#   rscgrp="interact"
   rscgrp="large"
fi

# Generate run.sh

cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# for FX10
#
################################################################################
#PJM --rsc-list "rscgrp=${rscgrp}"
#PJM --rsc-list "node=${TPROC}"
#PJM --rsc-list "elapse=12:00:00"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=16
export OMP_NUM_THREADS=16
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
