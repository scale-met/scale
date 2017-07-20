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
   SIN1_PP="#PJM --stgin  \"rank=* ${BINDIR}/${PPNAME}   %r:./\""

   CONFLIST=(`echo ${PPCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      SIN2_PP=`echo -e "${SIN2_PP}\n#PJM --stgin  \"rank=*         ./${CONFLIST[i]}   %r:./\""`
      RUN_PP=`echo -e "${RUN_PP}\n"${MPIEXEC} ${PROCLIST[i]} ./${PPNAME} ${CONFLIST[i]} "|| exit 1"`
   done
fi

if [ ! ${INITCONF} = "NONE" ]; then
   SIN1_INIT="#PJM --stgin  \"rank=* ${BINDIR}/${INITNAME} %r:./\""

   CONFLIST=(`echo ${INITCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      SIN2_INIT=`echo -e "${SIN2_INIT}\n#PJM --stgin  \"rank=*         ./${CONFLIST[i]}   %r:./\""`
      RUN_INIT=`echo -e "${RUN_INIT}\n"${MPIEXEC} ${PROCLIST[i]} ./${INITNAME} ${CONFLIST[i]} "|| exit 1"`
   done
fi

if [ ! ${RUNCONF} = "NONE" ]; then
   SIN1_MAIN="#PJM --stgin  \"rank=* ${BINDIR}/${BINNAME}  %r:./\""

   CONFLIST=(`echo ${RUNCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      SIN2_MAIN=`echo -e "${SIN2_MAIN}\n#PJM --stgin  \"rank=*         ./${CONFLIST[i]}   %r:./\""`
      RUN_MAIN=`echo -e "${RUN_MAIN}\n"fipp -C -Srange -Ihwm -d prof ${MPIEXEC} ${PROCLIST[i]} ./${BINNAME} ${CONFLIST[i]} "|| exit 1"`
   done
fi

if [ ! ${N2GCONF} = "NONE" ]; then
   SIN1_N2G="#PJM --stgin  \"rank=* ${BINDIR}/${N2GNAME}  %r:./\""

   CONFLIST=(`echo ${N2GCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      SIN2_N2G=`echo -e "${SIN2_N2G}\n#PJM --stgin  \"rank=*         ./${CONFLIST[i]}   %r:./\""`
      RUN_N2G=`echo -e "${RUN_N2G}\n"${MPIEXEC} ${PROCLIST[i]} ./${N2GNAME} ${CONFLIST[i]} "|| exit 1"`
   done
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
#PJM --stgout "rank=* %r:./prof/* ./prof/"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8

rm -rf ./prof

# run
${RUN_PP}
${RUN_INIT}
${RUN_MAIN}
${RUN_N2G}

################################################################################
EOF2
