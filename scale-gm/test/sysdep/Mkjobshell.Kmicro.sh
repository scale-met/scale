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

if [ ${TPROC} -gt 1152 ]; then
   rscgrp="invalid"
else
   rscgrp="micro"
fi

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# ------ For K computer
#
################################################################################
#PJM --rsc-list "rscgrp=${rscgrp}"
#PJM --rsc-list "node=${TPROC}"
#PJM --rsc-list "elapse=00:30:00"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
#export fu08bf=1
export XOS_MMM_L_ARENA_FREE=2

ln -svf ${TOPDIR}/bin/${BINNAME} .
ln -svf ${TOPDIR}/scale-gm/test/data/mnginfo/${MNGINFO} .
ln -svf ${TOPDIR}/scale-gm/test/data/grid/vgrid/${VGRID} .
EOF1

for f in $( ls ${TOPDIR}/scale-gm/test/data/grid/boundary/${dir2d} )
do
   echo "ln -svf ${TOPDIR}/scale-gm/test/data/grid/boundary/${dir2d}/${f} ." >> run.sh
done

cat << EOF2 >> run.sh

# run
${MPIEXEC} ./${BINNAME} nhm_driver.cnf || exit

################################################################################
EOF2


cat << EOFICO2LL1 > ico2ll.sh
#! /bin/bash -x
################################################################################
#
# ------ For K computer
#
################################################################################
#PJM --rsc-list "rscgrp=${rscgrp}"
#PJM --rsc-list "node=${TPROC}"
#PJM --rsc-list "elapse=00:30:00"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
#export fu08bf=1

ln -svf ${TOPDIR}/bin/gm_fio_ico2ll .
ln -svf ${TOPDIR}/scale-gm/test/data/mnginfo/${MNGINFO} .
ln -svf ${TOPDIR}/scale-gm/test/data/zaxis .
EOFICO2LL1

for f in $( ls ${TOPDIR}/scale-gm/test/data/grid/llmap/gl${GL}/rl${RL}/ )
do
   echo "ln -svf ${TOPDIR}/scale-gm/test/data/grid/llmap/gl${GL}/rl${RL}/${f} ." >> ico2ll.sh
done

cat << EOFICO2LL2 >> ico2ll.sh

# run
${MPIEXEC} ./gm_fio_ico2ll \
history \
glevel=${GLEV} \
rlevel=${RLEV} \
mnginfo="./${MNGINFO}" \
layerfile_dir="./zaxis" \
llmap_base="./llmap" \
-lon_swap \
-comm_smallchunk

################################################################################
EOFICO2LL2
