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
#PJM --rsc-list "elapse=01:00:00"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=16
export OMP_NUM_THREADS=16
#export fu08bf=1

rm -rf ./prof
mkdir -p ./prof

# run
                              ${MPIEXEC} ${BINDIR}/${INITNAME} ${INITCONF} || exit
fipp -C -Srange -Ihwm -d prof ${MPIEXEC} ${BINDIR}/${BINNAME}  ${RUNCONF}  || exit

################################################################################
EOF1
