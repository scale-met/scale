#! /bin/bash -x

GLEV=${1}
RLEV=${2}
NMPI=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

# System specific
MPIEXEC="mpiexec"

GL=`printf %02d ${GLEV}`
RL=`printf %02d ${RLEV}`
if   [ ${NMPI} -ge 10000 ]; then
	NP=`printf %05d ${NMPI}`
elif [ ${NMPI} -ge 1000 ]; then
	NP=`printf %04d ${NMPI}`
elif [ ${NMPI} -ge 100 ]; then
	NP=`printf %03d ${NMPI}`
else
	NP=`printf %02d ${NMPI}`
fi

dir2d=gl${GL}rl${RL}pe${NP}
dir3d=gl${GL}rl${RL}z${ZL}pe${NP}
res2d=GL${GL}RL${RL}
res3d=GL${GL}RL${RL}z${ZL}

MNGINFO=rl${RL}-prc${NP}.info

# for K(micro)
if [ ${NMPI} -gt 1152 ]; then
   rscgrp="invalid"
else
   rscgrp="micro"
fi

PROF1="fapp -C -Ihwm -Hevent=Cache        -d prof_cache -L 10"
PROF2="fapp -C -Ihwm -Hevent=Instructions -d prof_inst  -L 10"
PROF3="fapp -C -Ihwm -Hevent=MEM_access   -d prof_mem   -L 10"
PROF4="fapp -C -Ihwm -Hevent=Performance  -d prof_perf  -L 10"
PROF5="fapp -C -Ihwm -Hevent=Statistics   -d prof       -L 10"

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# for K micro
#
################################################################################
#PJM --rsc-list "rscgrp=${rscgrp}"
#PJM --rsc-list "node=${NMPI}"
#PJM --rsc-list "elapse=00:29:00"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
#export fu08bf=1
export XOS_MMM_L_ARENA_FREE=2

ln -sv ${TOPDIR}/bin/${BINNAME} .
ln -sv ${TOPDIR}/scale-gm/test/data/mnginfo/${MNGINFO} .
ln -sv ${TOPDIR}/scale-gm/test/data/grid/vgrid/${VGRID} .
EOF1

for f in $( ls ${TOPDIR}/scale-gm/test/data/grid/boundary/${dir2d} )
do
   echo "ln -sv ${TOPDIR}/scale-gm/test/data/grid/boundary/${dir2d}/${f} ." >> run.sh
done

cat << EOF2 >> run.sh
rm -rf ./prof*
mkdir -p ./prof_cache
mkdir -p ./prof_inst
mkdir -p ./prof_mem
mkdir -p ./prof_perf
mkdir -p ./prof

# run
${PROF1} ${MPIEXEC} ./${BINNAME} || exit
${PROF2} ${MPIEXEC} ./${BINNAME} || exit
${PROF3} ${MPIEXEC} ./${BINNAME} || exit
${PROF4} ${MPIEXEC} ./${BINNAME} || exit
${PROF5} ${MPIEXEC} ./${BINNAME} || exit

################################################################################
EOF2


cat << EOFICO2LL1 > ico2ll.sh
#! /bin/bash -x
################################################################################
#
# for K micro
#
################################################################################
#PJM --rsc-list "rscgrp=${rscgrp}"
#PJM --rsc-list "node=${NMPI}"
#PJM --rsc-list "elapse=00:29:00"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
#export fu08bf=1

ln -sv ${TOPDIR}/bin/gm_fio_ico2ll .
ln -sv ${TOPDIR}/scale-gm/test/data/mnginfo/${MNGINFO} .
ln -sv ${TOPDIR}/scale-gm/test/data/zaxis .
EOFICO2LL1

for f in $( ls ${TOPDIR}/scale-gm/test/data/grid/llmap/gl${GL}/rl${RL}/ )
do
   echo "ln -sv ${TOPDIR}/scale-gm/test/data/grid/llmap/gl${GL}/rl${RL}/${f} ." >> ico2ll.sh
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
