#! /bin/bash -x

GLEV=${1}
RLEV=${2}
TPROC=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# ------ FOR MacOSX & gfortran4.9 & OpenMPI1.7 -----
#
################################################################################
export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

ln -svf ${TOPDIR}/bin/${BINNAME} .

# run
./${BINNAME} mkmnginfo.cnf || exit

################################################################################
EOF1
