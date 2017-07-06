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
eval DATPARAM=(`echo ${9}  | tr -s '[' '"' | tr -s ']' '"'`)
eval DATDISTS=(`echo ${10} | tr -s '[' '"' | tr -s ']' '"'`)

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
#PJM --rsc-list "elapse=04:00:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
${SIN1_PP}
${SIN1_INIT}
${SIN1_MAIN}
${SIN2_PP}
${SIN2_INIT}
#PJM --stgin  "rank=* ./runA.conf  %r:./"
#PJM --stgin  "rank=* ./runB1.conf %r:./"
#PJM --stgin  "rank=* ./runB2.conf %r:./"
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

      pair=(${DATDISTS[$i]})

      if [ -f ${pair[0]}.pe000000.nc ]; then
         echo "#PJM --stgin  'rank=* ${pair[0]}.pe%06r.nc %r:./${pair[1]}.pe%06r.nc'" >> ./run.sh
      else
         echo "datafile does not found! : ${pair[0]}.pe000000.nc"
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
${MPIEXEC} ${BINDIR}/${BINNAME} runA.conf  || exit
${MPIEXEC} ${BINDIR}/${BINNAME} runB1.conf || exit
${MPIEXEC} ${BINDIR}/${BINNAME} runB2.conf || exit

################################################################################
EOF2