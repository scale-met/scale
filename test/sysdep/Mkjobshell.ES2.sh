#! /bin/bash -x

# Arguments
BINDIR=${1}
INITNAME=${2}
BINNAME=${3}
INITCONF=${4}
RUNCONF=${5}
TPROC=${6}

# System specific
MPIEXEC="mpirun -nnp ${TPROC}"

RUNDIR=`pwd`

# Generate run.sh

cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# for Earth Simulator 2 (L system)
#
################################################################################
#PBS -T mpisx
#PBS -q L
#PBS -b 1
#PBS -l elapstim_req=00:30:00
#PBS -l filecap_job=1gb

#PBS -v F_RECLUNIT=byte
#PBS -v F_ERRCNT=1
#PBS -v F_FTRACE=FMT1
#PBS -v F_PROGINF=DETAIL

#PBS -I "${BINDIR}/${INITNAME},ALL:./"
#PBS -I "${BINDIR}/${BINNAME},ALL:./"
#PBS -I "${RUNDIR}/${INITCONF},ALL:./"
#PBS -I "${RUNDIR}/${RUNCONF},ALL:./"

#PBS -O "${RUNDIR}/,0:./"

# run
${MPIEXEC} ./${INITNAME} ${INITCONF} || exit
${MPIEXEC} ./${BINNAME}  ${RUNCONF}  || exit

################################################################################
EOF1

cat << EOF1 > ./run_S.sh
#! /bin/bash -x
################################################################################
#
# for Earth Simulator 2 (S system)
#
################################################################################
#PBS -T mpisx
#PBS -q S
#PBS -l cpunum_job=1
#PBS -l cputim_job=00:30:00
#PBS -l memsz_job=60gb
#
#PBS -v F_RECLUNIT=byte
#PBS -v F_ERRCNT=1
#PBS -v F_FTRACE=FMT1
#PBS -v F_PROGINF=DETAIL

cd ${RUNDIR}
ln -v ${BINDIR}/${INITNAME} .
ln -v ${BINDIR}/${BINNAME}  .

# run
${MPIEXEC} ./${INITNAME} ${INITCONF} || exit
${MPIEXEC} ./${BINNAME}  ${RUNCONF}  || exit

################################################################################
EOF1
