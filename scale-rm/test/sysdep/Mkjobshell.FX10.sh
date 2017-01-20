#! /bin/bash -x#! /bin/bash -x

# Arguments
BINDIR=${1}
UTILDIR=${2}
PPNAME=${3}
INITNAME=${4}
BINNAME=${5}
N2GNAME=${6}
PPCONF=${7}
INITCONF=${8}
RUNCONF=${9}
N2GCONF=${10}
TPROC=${11}
DATDIR=${12}
DATPARAM=(`echo ${13} | tr -s ',' ' '`)
DATDISTS=(`echo ${14} | tr -s ',' ' '`)

# System specific
MPIEXEC="mpiexec"

if [ ! ${PPCONF} = "NONE" ]; then
  RUN_PP="${MPIEXEC} ${BINDIR}/${PPNAME} ${PPCONF} || exit"
fi

if [ ! ${INITCONF} = "NONE" ]; then
  RUN_INIT="${MPIEXEC} ${BINDIR}/${INITNAME} ${INITCONF} || exit"
fi

if [ ! ${RUNCONF} = "NONE" ]; then
  RUN_BIN="${MPIEXEC} ${BINDIR}/${BINNAME} ${RUNCONF} || exit"
fi

if [ ! ${N2GCONF} = "NONE" ]; then
  RUN_N2G="${MPIEXEC} ${UTILDIR}/${N2GNAME} ${N2GCONF} || exit"
fi

array=( `echo ${TPROC} | tr -s 'x' ' '`)
x=${array[0]}
y=${array[1]:-1}
let xy="${x} * ${y}"

# for AICS-FX10
if [ ${xy} -gt 96 ]; then
   rscgrp="huge"
elif [ ${xy} -gt 24 ]; then
   rscgrp="large"
else
   rscgrp="large"
fi





cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# ------ For FX10
#
################################################################################
#PJM --rsc-list "rscgrp=${rscgrp}"
#PJM --rsc-list "node=${TPROC}"
#PJM --rsc-list "elapse=12:00:00"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=16
export OMP_NUM_THREADS=16

EOF1

if [ ! ${DATPARAM[0]} = "" ]; then
   for f in ${DATPARAM[@]}
   do
         if [ -f ${DATDIR}/${f} ]; then
            echo "ln -svf ${DATDIR}/${f} ." >> ./run.sh
         elif [ -d ${DATDIR}/${f} ]; then
            echo "rm -f                  ./input" >> ./run.sh
            echo "ln -svf ${DATDIR}/${f} ./input" >> ./run.sh
         else
            echo "datafile does not found! : ${DATDIR}/${f}"
            exit 1
         fi
   done
fi

if [ ! ${DATDISTS[0]} = "" ]; then
   for prc in `seq 1 ${TPROC}`
   do
      let "prcm1 = ${prc} - 1"
      PE=`printf %06d ${prcm1}`
      for f in ${DATDISTS[@]}
      do
         if [ -f ${f}.pe${PE}.nc ]; then
            echo "ln -svf ${f}.pe${PE}.nc ." >> ./run.sh
         else
            echo "datafile does not found! : ${f}.pe${PE}.nc"
            exit 1
         fi
      done
   done
fi

cat << EOF2 >> ./run.sh

# run
${RUN_PP}
${RUN_INIT}
${RUN_BIN}
${RUN_N2G}

################################################################################
EOF2
