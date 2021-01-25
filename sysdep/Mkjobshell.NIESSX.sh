#!/bin/bash -x

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
MPIEXEC="mpirun"

PROCLIST=(`echo ${PROCS} | tr -s ',' ' '`)
TPROC=${PROCLIST[0]}
for n in ${PROCLIST[@]}
do
   (( n > TPROC )) && TPROC=${n}
done

NNODE=`expr \( $TPROC - 1 \) / 4 + 1`

if [ ! ${PPCONF} = "NONE" ]; then
   CONFLIST=(`echo ${PPCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      nn=`expr \( ${PROCLIST[i]} - 1 \) / 4 + 1`
      nnp=`expr ${PROCLIST[i]} / $NNODE`
      RUN_PP=`echo -e "${RUN_PP}\n"${MPIEXEC} -nn ${nn} -nnp ${nnp} /usr/lib/mpi/mpisep.sh ${BINDIR}/${PPNAME} ${CONFLIST[i]} "|| exit 1"`
   done
fi

if [ ! ${INITCONF} = "NONE" ]; then
   CONFLIST=(`echo ${INITCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      nn=`expr \( ${PROCLIST[i]} - 1 \) / 4 + 1`
      nnp=`expr ${PROCLIST[i]} / $NNODE`
      RUN_INIT=`echo -e "${RUN_INIT}\n"${MPIEXEC} -nn ${nn} -nnp ${nnp} /usr/lib/mpi/mpisep.sh ${BINDIR}/${INITNAME} ${CONFLIST[i]} "|| exit 1"`
   done
fi

if [ ! ${RUNCONF} = "NONE" ]; then
   CONFLIST=(`echo ${RUNCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      nn=`expr \( ${PROCLIST[i]} - 1 \) / 4 + 1`
      nnp=`expr ${PROCLIST[i]} / $NNODE`
      RUN_MAIN=`echo -e "${RUN_MAIN}\n"${MPIEXEC} -nn ${nn} -nnp ${nnp} /usr/lib/mpi/mpisep.sh ${BINDIR}/${BINNAME} ${CONFLIST[i]} "|| exit 1"`
   done
fi

if [ ! ${N2GCONF} = "NONE" ]; then
   CONFLIST=(`echo ${N2GCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      nn=`expr \( ${PROCLIST[i]} - 1 \) / 4 + 1`
      nnp=`expr ${PROCLIST[i]} / $NNODE`
      RUN_N2G=`echo -e "${RUN_N2G}\n"${MPIEXEC} -nn ${nn} -nnp ${nnp} /usr/lib/mpi/mpisep.sh ${BINDIR}/${N2GNAME} ${CONFLIST[i]} "|| exit 1"`
   done
fi

if [[ ${BINNAME} =~ ^scale-gm ]]; then
   nc=""
else
   nc=".nc"
fi



cat << EOF3 > ./run.sh
#!/bin/sh
################################################################################
#
# ------ For NIES SX system
#
################################################################################
#PBS -A rg108
#PBS -T mpisx
#PBS -q v_deb
#PBS -b ${NNODE}
#PBS -l elapstim_req=00:30:00
#PBS -l memsz_job=40gb
#PBS -v MPISEPSELECT=4

#PBS -v F_RECLUNIT=byte
#PBS -v F_ERRCNT=0
##PBS -v F_FTRACE=FMT1
#PBS -v F_PROGINF=DETAIL
#PBS -v MPIPROGINF=ALL_DETAIL
#PBS -v F_SETBUF=102400

cd `pwd`
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
         echo "ln -sf ${src} ./${dst}" >> ./run.sh
      elif [ -d ${src} ]; then
         echo "rm -f          ./${dst}" >> ./run.sh
         echo "ln -sf ${src} ./${dst}" >> ./run.sh
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
            echo "ln -sf ${src} ./${dst}" >> ./run.sh
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
${RUN_N2G}

################################################################################
EOF4
