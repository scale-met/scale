#! /bin/bash -x

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
  RUN_PP="${MPIEXEC} ./${PPNAME} ${PPCONF} || exit"
fi

if [ ! ${INITCONF} = "NONE" ]; then
  RUN_INIT="${MPIEXEC} ./${INITNAME} ${INITCONF} || exit"
fi

if [ ! ${RUNCONF} = "NONE" ]; then
  RUN_BIN="${MPIEXEC} ./${BINNAME} ${RUNCONF} || exit"
fi

if [ ! ${N2GCONF} = "NONE" ]; then
  RUN_N2G="${MPIEXEC} ${UTILDIR}/${N2GNAME} ${N2GCONF} || exit"
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
#PJM --stgin  "rank=* ${BINDIR}/${PPNAME}   %r:./"
#PJM --stgin  "rank=* ${BINDIR}/${INITNAME} %r:./"
#PJM --stgin  "rank=* ${BINDIR}/${BINNAME}  %r:./"
#PJM --stgin  "rank=* ${UTILDIR}/${N2GNAME} %r:./"
#PJM --stgin  "rank=*         ./${PPCONF}   %r:./"
#PJM --stgin  "rank=*         ./${INITCONF} %r:./"
#PJM --stgin  "rank=*         ./${RUNCONF}  %r:./"
#PJM --stgin  "rank=*         ./${N2GCONF}  %r:./"
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
${RUN_BIN}
${RUN_N2G}

################################################################################
EOF2