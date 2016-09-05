#! /bin/bash -x

GLEV=${1}
RLEV=${2}
TPROC=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

GL=`printf %02d ${GLEV}`
RL=`printf %02d ${RLEV}`
if   [ ${TPROC} -ge 10000 ]; then
	NP=`printf %05d ${TPROC}`
elif [ ${TPROC} -ge 1000 ]; then
	NP=`printf %04d ${TPROC}`
elif [ ${TPROC} -ge 100 ]; then
	NP=`printf %03d ${TPROC}`
else
	NP=`printf %02d ${TPROC}`
fi

dir2d=gl${GL}rl${RL}pe${NP}
dir3d=gl${GL}rl${RL}z${ZL}pe${NP}
res2d=GL${GL}RL${RL}
res3d=GL${GL}RL${RL}z${ZL}

MNGINFO=rl${RL}-prc${NP}.info

# System specific
MPIEXEC="mpirun -nnp 2"

NNODE=`expr $TPROC / 2`
RUNDIR=`pwd`

cat << EOFS1 > ./run.sh
#!/bin/sh
################################################################################
#
# for Earth Simulator 3 (S system)
#
################################################################################
#PBS -q S
#PBS -b ${NNODE}
#PBS -l elapstim_req=00:30:00
#PBS -l filecap_job=32gb

#PBS -v F_RECLUNIT=byte
#PBS -v F_ERRCNT=0
#PBS -v F_PROGINF=DETAIL
#PBS -v MPIPROGINF=ALL_DETAIL
#PBS -v F_SETBUF=102400

cd ${RUNDIR}

ln -s ${TOPDIR}/bin/${BINNAME} .
ln -s ${TOPDIR}/scale-gm/test/data/mnginfo/${MNGINFO} .
ln -s ${TOPDIR}/scale-gm/test/data/grid/vgrid/${VGRID} .
EOFS1

for f in $( ls ${TOPDIR}/scale-gm/test/data/grid/boundary/${dir2d} )
do
   echo "ln -s ${TOPDIR}/scale-gm/test/data/grid/boundary/${dir2d}/${f} ." >> run.sh
done

cat << EOFS2 >> run.sh
# run
${MPIEXEC} ./${BINNAME} nhm_driver.cnf || exit

################################################################################
EOFS2

cat << EOFL1 > run_L.sh
#!/bin/sh
################################################################################
#
# for Earth Simulator 3 (L system)
#
################################################################################
#PBS -q L
#PBS -b ${NNODE}
#PBS -l elapstim_req=00:30:00
#PBS -l filecap_job=32gb

#PBS -I "${TOPDIR}/bin/${BINNAME},ALL:./"
#PBS -I "${RUNDIR}/nhm_driver.cnf,ALL:./"
#PBS -I "${TOPDIR}/scale-gm/test/data/mnginfo/${MNGINFO},ALL:./"
#PBS -I "${TOPDIR}/scale-gm/test/data/grid/vgrid/${VGRID},ALL:./"
#PBS -I "${TOPDIR}/scale-gm/test/data/grid/boundary/${dir2d}/boundary_${res2d}.pe%06r,ALL:./"
#PBS -O "${RUNDIR}/,ALL:./"

# run
${MPIEXEC} ./${BINNAME} nhm_driver.cnf || exit

################################################################################
EOFL1

