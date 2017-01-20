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
MPIEXEC="mpiexec"

if [ ! ${PPCONF} = "NONE" ]; then
  SIN1_PP="#PJM --stgin  \"rank=* ${BINDIR}/${PPNAME}   %r:./\""
  SIN2_PP="#PJM --stgin  \"rank=*         ./${PPCONF}   %r:./\""
  RUN_PP="${MPIEXEC} ./${PPNAME} ${PPCONF} || exit"
fi

if [ ! ${INITCONF} = "NONE" ]; then
  SIN1_INIT="#PJM --stgin  \"rank=* ${BINDIR}/${INITNAME} %r:./\""
  SIN2_INIT="#PJM --stgin  \"rank=*         ./${INITCONF} %r:./\""
  RUN_INIT="${MPIEXEC} ./${INITNAME} ${INITCONF} || exit"
fi

if [ ! ${RUNCONF} = "NONE" ]; then
  SIN1_MAIN="#PJM --stgin  \"rank=* ${BINDIR}/${BINNAME}  %r:./\""
  SIN2_MAIN="#PJM --stgin  \"rank=*         ./${RUNCONF}  %r:./\""
  RUN_MAIN="${MPIEXEC} ./${BINNAME} ${RUNCONF} || exit"
fi

array=( `echo ${TPROC} | tr -s 'x' ' '`)
x=${array[0]}
y=${array[1]:-1}
let xy="${x} * ${y}"

if [ ${xy} -gt 36864 ]; then
   rscgrp="huge"
elif [ ${xy} -gt 384 ]; then
   rscgrp="large"
else
   rscgrp="small"
fi





cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# ------ For K computer
#
################################################################################
#PJM --rsc-list "rscgrp=${rscgrp}"
#PJM --rsc-list "node=${TPROC}"
#PJM --rsc-list "elapse=02:00:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
${SIN1_PP}
${SIN1_INIT}
${SIN1_MAIN}
${SIN2_PP}
${SIN2_INIT}
${SIN2_MAIN}
EOF1

if [ ! ${DATPARAM[0]} = "" ]; then
   for f in ${DATPARAM[@]}
   do
         if [ -f ${DATDIR}/${f} ]; then
            echo "#PJM --stgin  'rank=* ${DATDIR}/${f}   %r:./'"       >> ./run.sh
         elif [ -d ${DATDIR}/${f} ]; then
            echo "#PJM --stgin  'rank=* ${DATDIR}/${f}/* %r:./input/'" >> ./run.sh
         else
            echo "datafile does not found! : ${DATDIR}/${f}"
            exit 1
         fi
   done
fi

if [ ! ${DATDISTS[0]} = "" ]; then
   for f in ${DATDISTS[@]}
   do
      if [ -f ${f}.pe000000.nc ]; then
         echo "#PJM --stgin  'rank=* ${f}.pe%06r.nc %r:./'" >> ./run.sh
      else
         echo "datafile does not found! : ${f}.pe000000.nc"
         exit 1
      fi
   done
fi

cat << EOF2 >> ./run.sh
#PJM --stgout "rank=* %r:./* ./"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8

# run
${RUN_PP}
${RUN_INIT}
${RUN_MAIN}

################################################################################
EOF2
