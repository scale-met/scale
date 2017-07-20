#! /bin/bash -x

GLEV=${1}
RLEV=${2}
TPROC=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# ------ For SGI ICE X (Linux64 & gnu fortran&C & openmpi + Torque -----
#
################################################################################
#PBS -q s
#PBS -l nodes=1
#PBS -N mkmnginfo
#PBS -l walltime=00:10:00
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

# run
./${BINNAME} mkmnginfo.cnf || exit

################################################################################
EOF1
