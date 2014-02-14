#! /bin/bash -x

# Arguments
BINDIR=${1}
INITNAME=${2}
BINNAME=${3}
INITCONF=${4}
RUNCONF=${5}
TPROC=${6}

# System specific
MPIEXEC="impijob"

# Generate run.sh

cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# ------ FOR Linux64 & intel C&fortran & intel mpi & LSF-----
#
################################################################################
#BSUB -n ${TPROC}
#BSUB -q dl
#BSUB -a intelmpi
#BSUB -J ${BINNAME}
#BSUB -o STDOUT
#BSUB -e STDERR
export FORT_FMT_RECL=400
export OMP_NUM_THREADS=1

cd ${RUNDIR}

# run
${MPIEXEC} ${BINDIR}/${INITNAME} ${INITCONF} || exit
${MPIEXEC} ${BINDIR}/${BINNAME}  ${RUNCONF}  || exit

################################################################################
EOF1
