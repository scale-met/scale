#! /bin/bash -x

GLEV=${1}
RLEV=${2}
TPROC=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

# System specific
MPIEXEC="mpirun --mca btl openib,sm,self --bind-to core"

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

NNODE=`expr \( $TPROC - 1 \) / 20 + 1`
NPROC=`expr $TPROC / $NNODE`

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
# ------ For SGI ICE X (Linux64 & gnu fortran&C & openmpi + Torque -----
#
################################################################################
#PBS -q ${rscgrp}
#PBS -l nodes=${NNODE}:ppn=${NPROC}
#PBS -N mkrawgrid
#PBS -l walltime=00:30:00
#PBS -o STDOUT
#PBS -e STDERR
export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

source /etc/profile.d/modules.sh
module unload mpt/2.12
module unload intelcompiler
module unload intelmpi
module unload hdf5
module unload netcdf4/4.3.3.1-intel
module unload netcdf4/fortran-4.4.2-intel
module load gcc/4.7.2
module load openmpi/1.10.1-gcc
module load hdf5/1.8.16
module load netcdf4/4.3.3.1
module load netcdf4/fortran-4.4.2

cd \$PBS_O_WORKDIR

ln -svf ${TOPDIR}/bin/${BINNAME} .
EOF1

if   [ -f ${TOPDIR}/scale-gm/test/data/mnginfo/${MNGINFO} ]; then
   echo "mnginfo file is found in default database"
   echo "ln -svf ${TOPDIR}/scale-gm/test/data/mnginfo/${MNGINFO} ." >> run.sh
elif [ -f ../../mkmnginfo/rl${RL}pe${NP}/${MNGINFO} ]; then
   echo "mnginfo file is found in test directory"
   echo "ln -svf ../../mkmnginfo/rl${RL}pe${NP}/${MNGINFO} ." >> run.sh
else
   echo "mnginfo file is not found!"
   exit 1
fi

cat << EOF2 >> run.sh

# run
${MPIEXEC} ./${BINNAME} mkrawgrid.cnf || exit

################################################################################
EOF2
