#! /bin/bash -x

GLEV=${1}
RLEV=${2}
TPROC=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

# System specific
MPIEXEC="mpiexec"

GL=`printf %02d ${GLEV}`
RL=`printf %02d ${RLEV}`
if   [ ${TPROC} -ge 10000 ]; then
	NP=`printf %05d ${TPROC}`
elif [ ${TPROC} -ge 1000 ]; then
	NP=`printf %04d ${TPROC}`
elif [ ${TPROC} -ge 100 ]; then
	NP=`printf %03d ${TPROC}`
else
	NP=`printf %02d ${TPROC}`
fi

MNGINFO=rl${RL}-prc${NP}.info

# for K computer
if [ ${TPROC} -gt 36864 ]; then
   rscgrp="huge"
elif [ ${TPROC} -gt 384 ]; then
   rscgrp="large"
else
   rscgrp="small"
fi

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
#PJM --rsc-list "rscgrp=${rscgrp}"
#PJM --rsc-list "node=${TPROC}"
#PJM --rsc-list "elapse=00:30:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin  "rank=* ${TOPDIR}/bin/${BINNAME} %r:./"
#PJM --stgin  "rank=* ./mkrawgrid.cnf          %r:./"
EOF1

if   [ -f ${TOPDIR}/scale-gm/test/data/mnginfo/${MNGINFO} ]; then
   echo "mnginfo file is found in default database"
   echo "#PJM --stgin  \"rank=* ${TOPDIR}/scale-gm/test/data/mnginfo/${MNGINFO} %r:./\"" >> run.sh
elif [ -f ../../mkmnginfo/rl${RL}pe${NP}/${MNGINFO} ]; then
q   echo "mnginfo file is found in test directory"
   echo "#PJM --stgin  \"rank=* ../../mkmnginfo/rl${RL}pe${NP}/${MNGINFO} %r:./\"" >> run.sh
else
   echo "mnginfo file is not found!"
   exit 1
fi

cat << EOF2 >> run.sh
#PJM --stgout "rank=* %r:./*                      ./"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
export XOS_MMM_L_ARENA_FREE=2

# run
${MPIEXEC} ./${BINNAME} mkrawgrid.cnf || exit

################################################################################
EOF2
