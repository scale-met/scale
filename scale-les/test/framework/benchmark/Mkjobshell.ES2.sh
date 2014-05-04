#!/bin/bash -x

# Arguments
BINDIR=${1}
INITNAME=${2}
BINNAME=${3}
INITCONF=${4}
RUNCONF=${5}
TPROC=${6}

# System specific
MPIEXEC="mpirun -nnp ${TPROC} /usr/lib/mpi/mpisep.sh"

RUNDIR=`pwd`

# Generate run.sh

cat << EOF1 > ./run.sh
#!/bin/sh
################################################################################
#
# for Earth Simulator 2 (L system)
#
################################################################################
#PBS -T mpisx
#PBS -q L
#PBS -b 1
#PBS -l elapstim_req=01:00:00
#PBS -l filecap_job=100gb

#PBS -v F_RECLUNIT=byte
#PBS -v F_ERRCNT=0
##PBS -v F_FTRACE=FMT1
#PBS -v F_PROGINF=DETAIL
#PBS -v MPIPROGINF=ALL_DETAIL
#PBS -v F_SETBUF=102400
#PBS -v MPISEPSELECT=3

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

cat << EOF2 > ./run_S.sh
#!/bin/sh
################################################################################
#
# for Earth Simulator 2 (S system)
#
################################################################################
#PBS -T mpisx
#PBS -q S
#PBS -l cpunum_job=1
#PBS -l cputim_job=01:00:00
#PBS -l memsz_job=80gb

#PBS -v F_RECLUNIT=byte
#PBS -v F_ERRCNT=0
#PBS -v F_FTRACE=FMT1
#PBS -v F_PROGINF=DETAIL
#PBS -v MPIPROGINF=ALL_DETAIL
#PBS -v F_SETBUF=102400

cd ${RUNDIR}

# run
${MPIEXEC} ${BINDIR}/${INITNAME} ${INITCONF} || exit
${MPIEXEC} ${BINDIR}/${BINNAME}  ${RUNCONF}  || exit

################################################################################
EOF2
