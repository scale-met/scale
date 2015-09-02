#! /bin/bash -x

# Arguments
BINDIR=${1}
INITNAME=${2}
BINNAME=${3}
INITCONF=${4}
RUNCONF=${5}
TPROC=${6}
DATDIR=${7}
DATPARAM=(`echo ${8} | tr -s ',' ' '`)
DATDISTS=(`echo ${9} | tr -s ',' ' '`)

# System specific
MPIEXEC="mpiexec"

array=( `echo ${TPROC} | tr -s 'x' ' '`)
x=${array[0]}
y=${array[1]:-1}
let xy="${x} * ${y}"

if [ ${xy} -gt 1024 ]; then
   echo "Node usage is less than 1024. STOP"
   exit
fi

# PROF1="fapp -C -Ihwm -Hevent=Cache        -d prof_cache -L 10"
# PROF2="fapp -C -Ihwm -Hevent=Instructions -d prof_inst  -L 10"
# PROF3="fapp -C -Ihwm -Hevent=MEM_access   -d prof_mem   -L 10"
# PROF4="fapp -C -Ihwm -Hevent=Performance  -d prof_perf  -L 10"
PROF5="fapp -C -Ihwm -Hevent=Statistics   -d prof       -L 10"

# Generate run.sh

cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################
#PJM --rsc-list "rscgrp=micro"
#PJM --rsc-list "node=${TPROC}"
#PJM --rsc-list "elapse=00:05:00"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
#export fu08bf=1

# rm -rf prof_cache
# rm -rf prof_inst
# rm -rf prof_mem
# rm -rf prof_perf
rm -rf prof

EOF1

if [ ! ${DATPARAM[0]} = "" ]; then
   for f in ${DATPARAM[@]}
   do
         if [ -f ${DATDIR}/${f} ]; then
            echo "ln -svf ${DATDIR}/${f} ." >> ./run.sh
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
         if [ -f ${DATDIR}/${f}.pe${PE} ]; then
            echo "ln -svf ${DATDIR}/${f}.pe${PE} ." >> ./run.sh
         else
            echo "datafile does not found! : ${DATDIR}/${f}.pe${PE}"
            exit 1
         fi
      done
   done
fi

if [ ! ${PPNAME}   = "NONE" ]; then
   echo "${MPIEXEC} ${BINDIR}/${PPNAME}   ${PPCONF}   || exit 1" >> ./run.sh
fi

if [ ! ${INITNAME} = "NONE" ]; then
   echo "${MPIEXEC} ${BINDIR}/${INITNAME} ${INITCONF} || exit 1" >> ./run.sh
fi

if [ ! ${BINNAME} = "NONE" ]; then
#    echo "${PROF1} ${MPIEXEC} ${BINDIR}/${BINNAME}  ${RUNCONF}  || exit 1" >> ./run.sh
#    echo "${PROF2} ${MPIEXEC} ${BINDIR}/${BINNAME}  ${RUNCONF}  || exit 1" >> ./run.sh
#    echo "${PROF3} ${MPIEXEC} ${BINDIR}/${BINNAME}  ${RUNCONF}  || exit 1" >> ./run.sh
#    echo "${PROF4} ${MPIEXEC} ${BINDIR}/${BINNAME}  ${RUNCONF}  || exit 1" >> ./run.sh
   echo "${PROF5} ${MPIEXEC} ${BINDIR}/${BINNAME}  ${RUNCONF}  || exit 1" >> ./run.sh
fi
