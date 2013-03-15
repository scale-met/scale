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
y=${array[1]}
let xy="${x} * ${y}"

if [ ${xy} -gt 480 ]; then
   rscgrp="x-large"
elif [ ${xy} -gt 372 ]; then
   rscgrp="large"
elif [ ${xy} -gt 216 ]; then
   rscgrp="medium"
elif [ ${xy} -gt 12 ]; then
   rscgrp="small"
else
   rscgrp="short"
fi

# Generate run.sh

cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# for OAKLEAF-FX
#
################################################################################
#PJM --rsc-list "rscgrp=${rscgrp}"
#PJM --rsc-list "node=${TPROC}"
#PJM --rsc-list "elapse=00:30:00"
#PJM -j
#PJM -s
export PARALLEL=8
export OMP_NUM_THREADS=8
export fu08bf=10

ln -sv ${BINDIR}/${INITNAME} .
ln -sv ${BINDIR}/${BINNAME}  .

fprof="fipp -C -Srange -Ihwm -d prof"
rm -rf ./prof
mkdir -p ./prof

# run
          ${MPIEXEC} ./${INITNAME} ${INITCONF} || exit
\${fprof} ${MPIEXEC} ./${BINNAME}  ${RUNCONF}  || exit

################################################################################
EOF1
