#! /bin/bash -x

# Arguments
BINDIR=${1}
INITNAME=${2}
BINNAME=${3}
INITCONF=${4}
RUNCONF=${5}
TPROC=${6}

# System specific
MPIEXEC="openmpirun -np ${TPROC}"

# Generate run.sh

cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# for OAKLEAF-FX
#
################################################################################
#PJM --rsc-list "rscgrp=short"
#PJM --rsc-list "node=${TPROC}"
#PJM --rsc-list "elapse=00:30:00"
#PJM -j
#PJM -s
export PARALLEL=8
export OMP_NUM_THREADS=8
export LPG="/opt/FJSVxosmmm/sbin/lpgparm -s 32MB -d 32MB -h 32MB -t 32MB -p 32MB"
export FLIB_FASTOMP=TRUE
export THREAD_STACK_SIZE=8192
export XOS_MMM_L_ARENA_FREE=2

ln -sv ${BINDIR}/${INITNAME} .
ln -sv ${BINDIR}/${BINNAME}  .

export fprof="fipp -C -Ihwm -d prof"
rm -rf ./prof

# run
          ${MPIEXEC} ${INITNAME} ${INITCONF} || exit
\${fprof} ${MPIEXEC} ${BINNAME}  ${RUNCONF}  || exit

################################################################################
EOF1
