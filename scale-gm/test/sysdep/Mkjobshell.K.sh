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

dir2d=gl${GL}rl${RL}pe${NP}
dir3d=gl${GL}rl${RL}z${ZL}pe${NP}
res2d=GL${GL}RL${RL}
res3d=GL${GL}RL${RL}z${ZL}

MNGINFO=rl${RL}-prc${NP}.info

# for K computer
if [ ${TPROC} -gt 36864 ]; then
   rscgrp="huge"
elif [ ${TPROC} -gt 384 ]; then
   rscgrp="large"
else
   rscgrp="small"
fi

PROF1="fipp -C -Srange -Ihwm,nocall -d prof"
PROF2="fipp -C -Srange -Inohwm,call -d prof_call"

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
#PJM --stgin  "rank=* ${TOPDIR}/bin/${BINNAME}           %r:./"
#PJM --stgin  "rank=* ./nhm_driver.cnf                   %r:./"
#PJM --stgin  "rank=* ${TOPDIR}/scale-gm/test/data/mnginfo/${MNGINFO}  %r:./"
#PJM --stgin  "rank=* ${TOPDIR}/scale-gm/test/data/grid/vgrid/${VGRID} %r:./"
#PJM --stgin  "rank=* ${TOPDIR}/scale-gm/test/data/grid/boundary/${dir2d}/boundary_${res2d}.pe%06r %r:./"
#PJM --stgout "rank=* %r:./*           ./"
#PJM --stgout "rank=* %r:./prof/*      ./prof/"
#PJM --stgout "rank=* %r:./prof_call/* ./prof_call/"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
#export fu08bf=1
export XOS_MMM_L_ARENA_FREE=2

rm -rf ./prof
rm -rf ./prof_call
mkdir -p ./prof
mkdir -p ./prof_call

# run
${PROF1} ${MPIEXEC} ./${BINNAME} nhm_driver.cnf || exit
${PROF2} ${MPIEXEC} ./${BINNAME} nhm_driver.cnf || exit

################################################################################
EOF1


cat << EOFICO2LL1 > ico2ll.sh
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
#PJM --stgin  "rank=* ${TOPDIR}/bin/gm_fio_ico2ll      %r:./"
#PJM --stgin  "rank=* ${TOPDIR}/scale-gm/test/data/mnginfo/${MNGINFO} %r:./"
#PJM --stgin  "rank=* ${TOPDIR}/scale-gm/test/data/zaxis/*            %r:./"
#PJM --stgin  "rank=* ./history.pe%06r                  %r:./"
#PJM --stgin  "rank=* ${TOPDIR}/scale-gm/test/data/grid/llmap/gl${GL}/rl${RL}/llmap.* %r:./"
#PJM --stgout "rank=* %r:./*           ./"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
export fu08bf=1

# run
${MPIEXEC} ./gm_fio_ico2ll \
history \
glevel=${GLEV} \
rlevel=${RLEV} \
mnginfo="./${MNGINFO}" \
layerfile_dir="./." \
llmap_base="./llmap" \
-lon_swap \
-comm_smallchunk

################################################################################
EOFICO2LL1
