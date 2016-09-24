#!/bin/bash -x

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
MPIEXEC="mpirun -nnp ${TPROC} /usr/lib/mpi/mpisep.sh"

if [ ! ${PPCONF} = "NONE" ]; then
  RUN_PP="${MPIEXEC} ./${PPNAME} ${PPCONF} || exit"
fi

if [ ! ${INITCONF} = "NONE" ]; then
  RUN_INIT="${MPIEXEC} ./${INITNAME} ${INITCONF} || exit"
fi

if [ ! ${RUNCONF} = "NONE" ]; then
  RUN_BIN="${MPIEXEC} ./${BINNAME} ${RUNCONF} || exit"
fi

RUNDIR=`pwd`





cat << EOF1 > ./run.sh
#!/bin/sh
################################################################################
#
# ------ For Earth Simulator 3 (L system)
#
################################################################################
#PBS -T mpisx
#PBS -q L
#PBS -b 1
#PBS -l elapstim_req=01:00:00
#PBS -l filecap_job=100gb

#PBS -v F_RECLUNIT=byte
#PBS -v F_ERRCNT=0
#PBS -v F_PROGINF=DETAIL
#PBS -v MPIPROGINF=ALL_DETAIL
#PBS -v F_SETBUF=102400
#PBS -v MPISEPSELECT=3

#PBS -I "${BINDIR}/${PPNAME},ALL:./"
#PBS -I "${BINDIR}/${INITNAME},ALL:./"
#PBS -I "${BINDIR}/${BINNAME},ALL:./"
#PBS -I "${RUNDIR}/${PPCONF},ALL:./"
#PBS -I "${RUNDIR}/${INITCONF},ALL:./"
#PBS -I "${RUNDIR}/${RUNCONF},ALL:./"

#PBS -O "${RUNDIR}/,0:./"

# run
${RUN_PP}
${RUN_INIT}
${RUN_BIN}

################################################################################
EOF1





if [ ! ${PPCONF} = "NONE" ]; then
  RUN_PP="${MPIEXEC} ${BINDIR}/${PPNAME} ${PPCONF} || exit"
fi

if [ ! ${INITCONF} = "NONE" ]; then
  RUN_INIT="${MPIEXEC} ${BINDIR}/${INITNAME} ${INITCONF} || exit"
fi

if [ ! ${RUNCONF} = "NONE" ]; then
  RUN_BIN="${MPIEXEC} ${BINDIR}/${BINNAME} ${RUNCONF} || exit"
fi

cat << EOF2 > ./run_S.sh
#!/bin/sh
################################################################################
#
# ------ For Earth Simulator 3 (S system)
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
${RUN_PP}
${RUN_INIT}
${RUN_BIN}

################################################################################
EOF2
