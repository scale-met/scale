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
  RUN_PP="${MPIEXEC} ${BINDIR}/${PPNAME} ${PPCONF} || exit"
fi

if [ ! ${INITCONF} = "NONE" ]; then
  RUN_INIT="${MPIEXEC} ${BINDIR}/${INITNAME} ${INITCONF} || exit"
fi

if [ ! ${RUNCONF} = "NONE" ]; then
  RUN_BIN="fipp -C -Srange -Ihwm -d prof ${MPIEXEC} ${BINDIR}/${BINNAME} ${RUNCONF} || exit"
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
         echo "ln -svf ${src} ./${dst}" >> ./run.sh
      elif [ -d ${src} ]; then
         echo "rm -f          ./${dst}" >> ./run.sh
         echo "ln -svf ${src} ./${dst}" >> ./run.sh
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

      for np in `seq 1 ${triple[0]}`
      do
         let "ip = ${np} - 1"
         PE=`printf %06d ${ip}`

         src=${triple[1]}.pe${PE}.nc
         dst=${triple[2]}.pe${PE}.nc

         if [ -f ${src} ]; then
            echo "ln -svf ${src} ./${dst}" >> ./run.sh
         else
            echo "datafile does not found! : ${src}"
            exit 1
         fi
      done
   done
fi

cat << EOF2 >> ./run.sh

rm -rf ./prof
mkdir -p ./prof

# run
${RUN_PP}
${RUN_INIT}
${RUN_BIN}

################################################################################
EOF2
