#! /bin/bash -x

# Arguments
BINDIR=${1}
PPNAME=${2}
INITNAME=${3}
BINNAME=${4}
PPCONF=${5}
INITCONF=${6}
RUNCONF=${7}
TPROC=${8}
DATDIR=${9}
DATPARAM=(`echo ${10} | tr -s ',' ' '`)
DATDISTS=(`echo ${11} | tr -s ',' ' '`)

# System specific
MPIEXEC="mpiexec_mpt -n ${TPROC}"

# Generate run.sh

cat << EOF1 > ./run.sh
#! /bin/bash
#PBS -l nodes=1:ppn=4
##PBS -l select=${TPROC}:ncpus=1
################################################################################
#
# ------ FOR Linux64 & intel C&fortran & intel SGI mpt -----
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


#rm -f             ./input
#ln -svf ${DATDIR} ./input

# run
EOF1

if [ ! ${PPNAME}   = "NONE" ]; then
   echo "${MPIEXEC} ${BINDIR}/${PPNAME}   ${PPCONF}   || exit 1" >> ./run.sh
fi

if [ ! ${INITNAME} = "NONE" ]; then
   echo "${MPIEXEC} ${BINDIR}/${INITNAME} ${INITCONF} || exit 1" >> ./run.sh
fi

if [ ! ${BINNAME} = "NONE" ]; then
   echo "${MPIEXEC} ${BINDIR}/${BINNAME}  ${RUNCONF}  || exit 1" >> ./run.sh
fi
