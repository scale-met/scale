#! /bin/bash -x

# Arguments
BINDIR=${1}
PPNAME=${2}
INITNAME=${3}
BINNAME=${4}
N2GNAME=${5}
PPCONF=${6}
INITCONF=${7}
RUNCONF=${8}
N2GCONF=${9}
PROCS=${10}
eval DATPARAM=(`echo ${11} | tr -s '[' '"' | tr -s ']' '"'`)
eval DATDISTS=(`echo ${12} | tr -s '[' '"' | tr -s ']' '"'`)

# System specific
MPIEXEC="mpiexec -np"

PROCLIST=(`echo ${PROCS} | tr -s ',' ' '`)
TPROC=${PROCLIST[0]}
for n in ${PROCLIST[@]}
do
   (( n > TPROC )) && TPROC=${n}
done

if [ ! ${PPCONF} = "NONE" ]; then
   CONFLIST=(`echo ${PPCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      RUN_PP=`echo -e "${RUN_PP}\n"${MPIEXEC} ${PROCLIST[i]} ${BINDIR}/${PPNAME} ${CONFLIST[i]} "|| exit 1"`
   done
fi

if [ ! ${INITCONF} = "NONE" ]; then
   CONFLIST=(`echo ${INITCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      RUN_INIT=`echo -e "${RUN_INIT}\n"${MPIEXEC} ${PROCLIST[i]} ${BINDIR}/${INITNAME} ${CONFLIST[i]} "|| exit 1"`
   done
fi

if [ ! ${RUNCONF} = "NONE" ]; then
   CONFLIST=(`echo ${RUNCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      RUN_MAIN=`echo -e "${RUN_MAIN}\n"${MPIEXEC} ${PROCLIST[i]} ${BINDIR}/${BINNAME} ${CONFLIST[i]} "|| exit 1"`
   done
fi

if [ ! ${N2GCONF} = "NONE" ]; then
   CONFLIST=(`echo ${N2GCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      RUN_N2G=`echo -e "${RUN_N2G}\n"${MPIEXEC} ${PROCLIST[i]} ${BINDIR}/${N2GNAME} ${CONFLIST[i]} "|| exit 1"`
   done
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

if [ "${BINNAME}" = "scale-gm" ]; then
   nc=""
else
   nc=".nc"
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

         src=${triple[1]}.pe${PE}${nc}
         dst=${triple[2]}.pe${PE}${nc}

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

# run
${RUN_PP}
${RUN_INIT}
${RUN_MAIN}
${RUN_N2G}

################################################################################
EOF2
