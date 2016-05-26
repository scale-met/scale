#! /bin/bash -x

# Arguments
BINDIR=${1}
INITNAME=${2}
BINNAME=${3}
INITCONF=${4}
RUNCONF=${5}
TPROC=${6}
DATDIR=${7}
DATPARAM=(`echo ${8} | tr -s ',' ' '`)
DATDISTS=(`echo ${9} | tr -s ',' ' '`)

# System specific
MPIEXEC="mpirun --mca btl openib,self -np ${TPROC}"

NNODE=`expr \( $TPROC - 1 \) / 8 + 1`
NPROC=`expr $TPROC / $NNODE`

if [ ${NNODE} -gt 16 ]; then
   rscgrp="l"
elif [ ${TPROC} -gt 3 ]; then
   rscgrp="m"
else
   rscgrp="s"
fi

# Generate run.sh

cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# ------ FOR SGI ICE X (Linux64 & gnu C&fortran & openmpi + Torque -----
#
################################################################################
#PBS -q ${rscgrp}
#PBS -l nodes=${NNODE}:ppn=${NPROC}
#PBS -l walltime=1:00:00
#PBS -N SCALE
#PBS -o OUT.log
#PBS -e ERR.log
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

EOF1

if [ ! ${DATPARAM[0]} = "" ]; then
   for f in ${DATPARAM[@]}
   do
         if [ -f ${DATDIR}/${f} ]; then
            echo "ln -svf ${DATDIR}/${f} ." >> ./run.sh
         else
            echo "datafile does not found! : ${DATDIR}/${f}"
            exit 1
         fi
   done
fi

if [ ! ${DATDISTS[0]} = "" ]; then
   for prc in `seq 1 ${TPROC}`
   do
      let "prcm1 = ${prc} - 1"
      PE=`printf %06d ${prcm1}`
      for f in ${DATDISTS[@]}
      do
         if [ -f ${f}.pe${PE}.nc ]; then
            echo "ln -svf ${f}.pe${PE}.nc ." >> ./run.sh
         else
            echo "datafile does not found! : ${f}.pe${PE}.nc"
            exit 1
         fi
      done
   done
fi


if [ ! ${INITNAME} = "NONE" ]; then
  RUN_INIT="${MPIEXEC} ${BINDIR}/${INITNAME} ${INITCONF} || exit"
fi

if [ ! ${BINNAME} = "NONE" ]; then
  RUN_BIN="${MPIEXEC} ${BINDIR}/${BINNAME} ${RUNCONF} || exit"
fi

cat << EOF2 >> ./run.sh

# run
${RUN_INIT}
${RUN_BIN}

################################################################################
EOF2
