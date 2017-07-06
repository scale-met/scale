#!/bin/bash -x

# Arguments
BINDIR=${1}
PPNAME=${2}
INITNAME=${3}
BINNAME=${4}
PPCONF=${5}
INITCONF=${6}
RUNCONF=${7}
PROCS=${8}
eval DATPARAM=(`echo ${9}  | tr -s '[' '"' | tr -s ']' '"'`)
eval DATDISTS=(`echo ${10} | tr -s '[' '"' | tr -s ']' '"'`)

# System specific
MPIEXEC="mpirun -nnp"

PROCLIST=(`echo ${PROCS} | tr -s ',' ' '`)
TPROC=${PROCLIST[0]}
for n in ${PROCLIST[@]}
do
   (( n > TPROC )) && TPROC=${n}
done

RUNDIR=`pwd`

if [ ! ${PPCONF} = "NONE" ]; then
  SIN1_PP="#PBS -I \"${BINDIR}/${PPNAME},ALL:./\""
  SIN2_PP="#PBS -I \"${RUNDIR}/${PPCONF},ALL:./\""
  RUN_PP="${MPIEXEC} ./${PPNAME} ${PPCONF} || exit"
fi

if [ ! ${INITCONF} = "NONE" ]; then
  SIN1_INIT="#PBS -I \"${BINDIR}/${INITNAME},ALL:./\""
  SIN2_INIT="#PBS -I \"${RUNDIR}/${INITCONF},ALL:./\""
  RUN_INIT="${MPIEXEC} ./${INITNAME} ${INITCONF} || exit"
fi

if [ ! ${RUNCONF} = "NONE" ]; then
  SIN1_MAIN="#PBS -I \"${BINDIR}/${BINNAME},ALL:./\""
  SIN2_MAIN="#PBS -I \"${RUNDIR}/${RUNCONF},ALL:./\""
  RUN_MAIN="${MPIEXEC} ./${BINNAME} ${RUNCONF} || exit"
fi





cat << EOF1 > ./run_L.sh
#!/bin/sh
################################################################################
#
# ------ For Earth Simulator 3 (L system)
#
################################################################################
#PBS -T mpisx
#PBS -q L
#PBS -b 1
#PBS -l elapstim_req=01:00:00
#PBS -l filecap_job=100gb

#PBS -v F_RECLUNIT=byte
#PBS -v F_ERRCNT=0
#PBS -v F_PROGINF=DETAIL
#PBS -v MPIPROGINF=ALL_DETAIL
#PBS -v F_SETBUF=102400
#PBS -v MPISEPSELECT=3
${SIN1_PP}
${SIN1_INIT}
${SIN1_MAIN}
${SIN2_PP}
${SIN2_INIT}
${SIN2_MAIN}
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
         echo "#PBS -I \"${src},ALL:./\"" >> ./run_L.sh
      elif [ -d ${src} ]; then
         echo "#PBS -I \"${src}/*,ALL:./\"" >> ./run_L.sh
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
         echo "#PBS -I \"${{triple[1]}.pe%06r.nc,ALL:./${triple[2]}.pe%06r.nc\"" >> ./run_L.sh
      else
         echo "datafile does not found! : ${{triple[1]}.pe000000.nc"
         exit 1
      fi
   done
fi

cat << EOF2 >> ./run_L.sh
#PBS -O "${RUNDIR}/,0:./"

# run
${RUN_PP}
${RUN_INIT}
${RUN_MAIN}

################################################################################
EOF2





if [ ! ${PPCONF} = "NONE" ]; then
   CONFLIST=(`echo ${PPCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      RUN_PP=`echo -e "${RUN_PP}\n"${MPIEXEC} ${PROCLIST[i]} ${BINDIR}/${PPNAME} ${CONFLIST[i]} || exit`
   done
fi

if [ ! ${INITCONF} = "NONE" ]; then
   CONFLIST=(`echo ${INITCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      RUN_INIT=`echo -e "${RUN_INIT}\n"${MPIEXEC} ${PROCLIST[i]} ${BINDIR}/${INITNAME} ${CONFLIST[i]} || exit`
   done
fi

if [ ! ${RUNCONF} = "NONE" ]; then
   CONFLIST=(`echo ${RUNCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      RUN_MAIN=`echo -e "${RUN_MAIN}\n"${MPIEXEC} ${PROCLIST[i]} ${BINDIR}/${BINNAME} ${CONFLIST[i]} || exit`
   done
fi

cat << EOF3 > ./run.sh
#!/bin/sh
################################################################################
#
# ------ For Earth Simulator 3 (S system)
#
################################################################################
#PBS -T mpisx
#PBS -q S
#PBS -l cpunum_job=1
#PBS -l cputim_job=01:00:00
#PBS -l memsz_job=80gb

#PBS -v F_RECLUNIT=byte
#PBS -v F_ERRCNT=0
##PBS -v F_FTRACE=FMT1
#PBS -v F_PROGINF=DETAIL
#PBS -v MPIPROGINF=ALL_DETAIL
#PBS -v F_SETBUF=102400

cd ${RUNDIR}
EOF3

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

cat << EOF4 >> ./run.sh

# run
${RUN_PP}
${RUN_INIT}
${RUN_MAIN}

################################################################################
EOF4
