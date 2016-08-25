#! /bin/bash -x

GLEV=${1}
RLEV=${2}
NMPI=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

# System specific
MPIEXEC="mpirun --mca btl openib,self"

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

NNODE=`expr \( $NMPI - 1 \) / 10 + 1`
NPROC=`expr $NMPI / $NNODE`

if [ ${NNODE} -gt 16 ]; then
   rscgrp="l"
elif [ ${NNODE} -gt 3 ]; then
   rscgrp="m"
else
   rscgrp="s"
fi

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# ------ For Linux64 & gnu fortran&C & OpenMPI & torque -----
#
################################################################################
#PBS -q ${rscgrp}
#PBS -l nodes=${NNODE}:ppn=${NPROC}
#PBS -N ${res3d}
#PBS -l walltime=04:00:00
#PBS -o STDOUT
#PBS -e STDERR
export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

cd \$PBS_O_WORKDIR

ln -sv ${TOPDIR}/bin/${BINNAME} .
ln -sv ${TOPDIR}/scale-gm/test/data/mnginfo/${MNGINFO} .
ln -sv ${TOPDIR}/scale-gm/test/data/grid/vgrid/${VGRID} .
EOF1

for f in $( ls ${TOPDIR}/scale-gm/test/data/grid/boundary/${dir2d} )
do
   echo "ln -sv ${TOPDIR}/scale-gm/test/data/grid/boundary/${dir2d}/${f} ." >> run.sh
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
# ------ For Linux64 & gnu fortran&C & OpenMPI & torque -----
#
################################################################################
#PBS -q ${rscgrp}
#PBS -l nodes=${NNODE}:ppn=${NPROC}
#PBS -N ico2ll_${res3d}
#PBS -l walltime=10:00:00
#PBS -o STDOUT
#PBS -e STDERR
export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

cd \$PBS_O_WORKDIR

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
