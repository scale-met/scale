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
MPIEXEC="mpiexec_mpt -n ${TPROC}"

# Generate run.sh

cat << EOF1 > ./run.sh
#! /bin/bash
#PBS -l nodes=1:ppn=4
##PBS -l select=${TPROC}:ncpus=1
################################################################################
#
# ------ FOR Linux64 & intel C&fortran & SGI mpt -----
#
################################################################################
export FORT_FMT_RECL=400

. /usr/share/modules/init/sh  ### this is required to use "module" command
module unload pgi-12.10
module load intel-2013.1.117
module load mpt/2.05
module load common_intel/hdf5/1.8.14
module load common_intel/netcdf/4.3.2
module load common_intel/netcdf-fortran/4.4.1

ulimit -s unlimited

nodelist=''
for inode in \`cat \$PBS_NODEFILE | sort | uniq\`; do
  if [ -z "\$nodelist" ]; then
    nodelist="\$inode"
  else
    nodelist="\$nodelist,\$inode"
  fi
done

export MPI_UNIVERSE="\$nodelist 16"
echo "MPI_UNIVERSE='\$MPI_UNIVERSE'"

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

if [ "${INITCONF}" != 'NONE' ]; then
  echo "time ${MPIEXEC} ${BINDIR}/${INITNAME} ${INITCONF} || exit" >> ./run.sh
fi
if [ "${RUNCONF}" != 'NONE' ]; then
  echo "time ${MPIEXEC} ${BINDIR}/${BINNAME}  ${RUNCONF}  || exit" >> ./run.sh
fi
