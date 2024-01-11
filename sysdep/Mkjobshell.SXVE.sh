#! /bin/bash -x

# Arguments
MPIEXEC=${1}
BINDIR=${2}
PPNAME=${3}
INITNAME=${4}
BINNAME=${5}
N2GNAME=${6}
PPCONF=${7}
INITCONF=${8}
RUNCONF=${9}
N2GCONF=${10}
PROCS=${11}
eval DATPARAM=(`echo ${12} | tr -s '[' '"' | tr -s ']' '"'`)
eval DATDISTS=(`echo ${13} | tr -s '[' '"' | tr -s ']' '"'`)

PROCLIST=(`echo ${PROCS} | tr -s ',' ' '`)
TPROC=${PROCLIST[0]}
for n in ${PROCLIST[@]}
do
   (( n > TPROC )) && TPROC=${n}
done

NVH=`expr \( $TPROC - 1 \) / 64   + 1`
NVE=`expr \( $TPROC - 1 \) / 8    + 1`
NVE=`expr \( $NVE   - 1 \) / $NVH + 1`

if [ ! ${PPCONF} = "NONE" ]; then
   CONFLIST=(`echo ${PPCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      nn=`expr \( ${PROCLIST[i]} - 1 \) / 8 + 1`
      nnp=`expr ${PROCLIST[i]} / ${nn}`
      RUN_PP=`echo -e "${RUN_PP}\n"${MPIEXEC} -nn ${nn} -nnp ${nnp} ${BINDIR}/${PPNAME} ${CONFLIST[i]} "|| exit 1"`
   done
fi

if [ ! ${INITCONF} = "NONE" ]; then
   CONFLIST=(`echo ${INITCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      nn=`expr \( ${PROCLIST[i]} - 1 \) / 8 + 1`
      nnp=`expr ${PROCLIST[i]} / ${nn}`
      RUN_INIT=`echo -e "${RUN_INIT}\n"${MPIEXEC} -nn ${nn} -nnp ${nnp} ${BINDIR}/${INITNAME} ${CONFLIST[i]} "|| exit 1"`
   done
fi

if [ ! ${RUNCONF} = "NONE" ]; then
   CONFLIST=(`echo ${RUNCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      nn=`expr \( ${PROCLIST[i]} - 1 \) / 8 + 1`
      nnp=`expr ${PROCLIST[i]} / ${nn}`
      RUN_MAIN=`echo -e "${RUN_MAIN}\n"${MPIEXEC} -nn ${nn} -nnp ${nnp} ${BINDIR}/${BINNAME} ${CONFLIST[i]} "|| exit 1"`
   done
fi

if [ ! ${N2GCONF} = "NONE" ]; then
   CONFLIST=(`echo ${N2GCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      nn=`expr \( ${PROCLIST[i]} - 1 \) / 8 + 1`
      nnp=`expr ${PROCLIST[i]} / ${nn}`
      RUN_N2G=`echo -e "${RUN_N2G}\n"${MPIEXEC} -nn ${nn} -nnp ${nnp} ${BINDIR}/${N2GNAME} ${CONFLIST[i]} "|| exit 1"`
   done
fi

if [[ ${BINNAME} =~ ^scale-gm ]]; then
   nc=""
else
   nc=".nc"
fi



cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# ------ For SX-Aurora TSUBASA
#
################################################################################
#PBS -q v_normal
#PBS -T necmpi
#PBS -l elapstim_req=1:00:00
#PBS -v OMP_NUM_THREADS=1
##PBS -v VE_PROGINF=YES
#PBS -v VE_FORT_UFMTENDIAN=ALL
#PBS -b ${NVH}
#PBS --venum-lhost=${NVE}

module load ncc/5.1.0
module load nfort/5.1.0
module load nec-mpi/latest
module load nlc/latest
module load netCDF-f_ve/4.5.4

export OMP_NUM_THREADS=1
cd \$PBS_O_WORKDIR

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
