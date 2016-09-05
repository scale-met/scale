#! /bin/bash -x

GLEV=${1}
RLEV=${2}
TPROC=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

# System specific
MPIEXEC="srun"

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

NNODE=`expr \( $TPROC - 1 \) / 20 + 1`
NPROC=`expr $TPROC / $NNODE`

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# for DKRZ Mistral with SLURM
#
################################################################################
#SBATCH --partition=compute
#SBATCH --job-name=NICAMDC
#SBATCH --nodes=${NNODE}
#SBATCH --ntasks-per-node=${NPROC}
#SBATCH --time=00:30:00
set -e

export OMP_NUM_THREADS=1

ln -svf ${TOPDIR}/bin/${BINNAME} .
ln -svf ${TOPDIR}/data/mnginfo/${MNGINFO} .
ln -svf ${TOPDIR}/data/grid/vgrid/${VGRID} .
EOF1

for f in $( ls ${TOPDIR}/data/grid/boundary/${dir2d} )
do
   echo "ln -svf ${TOPDIR}/data/grid/boundary/${dir2d}/${f} ." >> run.sh
done

cat << EOF2 >> run.sh

${MPIEXEC} ./${BINNAME} || exit 1

################################################################################
EOF2


cat << EOFICO2LL1 > ico2ll.sh
#! /bin/bash -x
################################################################################
#
# for DKRZ Mistral with SLURM
#
################################################################################
#SBATCH --partition=compute
#SBATCH --job-name=ICO2LL
#SBATCH --nodes=${NNODE}
#SBATCH --ntasks-per-node=${NPROC}
#SBATCH --time=00:30:00
set -e

export OMP_NUM_THREADS=1

ln -svf ${TOPDIR}/bin/fio_ico2ll_mpi .
ln -svf ${TOPDIR}/data/mnginfo/${MNGINFO} .
ln -svf ${TOPDIR}/data/zaxis .
EOFICO2LL1

for f in $( ls ${TOPDIR}/data/grid/llmap/gl${GL}/rl${RL}/ )
do
   echo "ln -svf ${TOPDIR}/data/grid/llmap/gl${GL}/rl${RL}/${f} ." >> ico2ll.sh
done

cat << EOFICO2LL2 >> ico2ll.sh

# run
${MPIEXEC} ./fio_ico2ll_mpi \
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
