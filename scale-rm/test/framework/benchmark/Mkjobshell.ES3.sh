#!/bin/bash -x

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
MPIEXEC="mpirun -nnp ${TPROC} /usr/lib/mpi/mpisep.sh"

if [ ! ${INITNAME} = "NONE" ]; then
  RUN_INIT="${MPIEXEC} ./${INITNAME} ${INITCONF} || exit"
fi

if [ ! ${BINNAME} = "NONE" ]; then
  RUN_BIN="${MPIEXEC} ./${BINNAME} ${RUNCONF} || exit"
fi





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

#PBS -I "${BINDIR}/${INITNAME},ALL:./"
#PBS -I "${BINDIR}/${BINNAME},ALL:./"
#PBS -I "${RUNDIR}/${INITCONF},ALL:./"
#PBS -I "${RUNDIR}/${RUNCONF},ALL:./"

#PBS -O "${RUNDIR}/,0:./"

# run
${RUN_INIT}
${RUN_BIN}

################################################################################
EOF1





RUNDIR=`pwd`

if [ ! ${INITNAME} = "NONE" ]; then
  RUN_INIT="${MPIEXEC} ${BINDIR}/${INITNAME} ${INITCONF} || exit"
fi

if [ ! ${BINNAME} = "NONE" ]; then
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
${RUN_INIT}
${RUN_BIN}

################################################################################
EOF2
