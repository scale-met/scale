#! /bin/bash -x

# Arguments
BINDIR=${1}
NET2GCONF=${2}
TPROC=${3}

# System specific
MPIEXEC="mpirun -np ${TPROC}"

# Generate run.sh

cat << EOF1 > ./net2g.sh
#! /bin/bash -x
################################################################################
#
# ------ For Linux64 & GNU fortran&C & openmpi -----
#
################################################################################
export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

# run
${MPIEXEC} ${BINDIR}/net2g ${NET2GCONF} || exit 1

EOF1
