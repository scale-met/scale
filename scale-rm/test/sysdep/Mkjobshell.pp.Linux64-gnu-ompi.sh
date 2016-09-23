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
MPIEXEC="mpirun -np ${TPROC}"

# Generate run.sh

cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# ------ For Linux64 & GNU fortran&C & openmpi -----
#
################################################################################
export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

rm -f             ./input
ln -svf ${DATDIR} ./input

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
