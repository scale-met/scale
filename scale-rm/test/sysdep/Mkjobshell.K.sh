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
eval DATPARAM=(`echo ${12} | tr -s '[' '"' | tr -s ']' '"'`)
eval DATDISTS=(`echo ${13} | tr -s '[' '"' | tr -s ']' '"'`)

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

if [ ! ${N2GCONF} = "NONE" ]; then
  SIN1_N2G="#PJM --stgin  \"rank=* ${UTILDIR}/${N2GNAME} %r:./\""
  SIN2_N2G="#PJM --stgin  \"rank=*         ./${N2GCONF}  %r:./\""
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
${SIN1_PP}
${SIN1_INIT}
${SIN1_MAIN}
${SIN1_N2G}
${SIN2_PP}
${SIN2_INIT}
${SIN2_MAIN}
${SIN2_N2G}
EOF1

# link to file or directory
ndata=${#DATPARAM[@]}

if [ ${ndata} -gt 0 ]; then
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"

      pair=(${DATPARAM[$i]})

      src=${pair[0]}
      dst=${pair[1]}
      if [ "${dst}" = "" ]; then
         dst=${pair[0]}
      fi

      if [ -f ${src} ]; then
         echo "#PJM --stgin  'rank=* ${src}   %r:./${dst} '" >> ./run.sh
      elif [ -d ${src} ]; then
         echo "#PJM --stgin  'rank=* ${src}/* %r:./${dst}/'" >> ./run.sh
      else
         echo "datafile does not found! : ${src}"
         exit 1
      fi
   done
fi

# link to distributed file
ndata=${#DATDISTS[@]}

if [ ${ndata} -gt 0 ]; then
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"

      triple=(${DATDISTS[$i]})

      if [ -f ${triple[1]}.pe000000.nc ]; then
         echo "#PJM --stgin  'rank=* ${triple[1]}.pe%06r.nc %r:./${triple[2]}.pe%06r.nc'" >> ./run.sh
      else
         echo "datafile does not found! : ${triple[1]}.pe000000.nc"
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
${RUN_N2G}

################################################################################
EOF2
