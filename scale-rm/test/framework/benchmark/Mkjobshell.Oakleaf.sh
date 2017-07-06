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

# for Oakleaf-FX
# if [ ${xy} -gt 480 ]; then
#    rscgrp="x-large"
# elif [ ${xy} -gt 372 ]; then
#    rscgrp="large"
# elif [ ${xy} -gt 216 ]; then
#    rscgrp="medium"
# elif [ ${xy} -gt 12 ]; then
#    rscgrp="small"
# else
#    rscgrp="short"
# fi
rscgrp="debug"





cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# ------ For Oakleaf-FX -----
#
################################################################################
#PJM --rsc-list "rscgrp=${rscgrp}"
#PJM --rsc-list "node=${TPROC}"
#PJM --rsc-list "elapse=00:20:00"
#PJM -j
#PJM -s
#
module load netCDF
module load netCDF-fortran
module load HDF5/1.8.9
module list
#
export PARALLEL=8
export OMP_NUM_THREADS=8

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
