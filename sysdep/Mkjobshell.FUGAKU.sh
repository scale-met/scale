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

FILES_LLIO=""
if [ ! ${PPCONF} = "NONE" ]; then
   CONFLIST=(`echo ${PPCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   FILES_LLIO=`echo -e ${FILES_LLIO} ${BINDIR}/${PPNAME}`
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      RUN_PP=`echo -e "${RUN_PP}\n"${MPIEXEC} ${PROCLIST[i]} ${BINDIR}/${PPNAME} ${CONFLIST[i]} "|| exit 1"`
      FILES_LLIO=`echo -e ${FILES_LLIO} ${CONFLIST[i]}`
   done
fi

if [ ! ${INITCONF} = "NONE" ]; then
   CONFLIST=(`echo ${INITCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   FILES_LLIO=`echo -e ${FILES_LLIO} ${BINDIR}/${INITNAME}`
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      RUN_INIT=`echo -e "${RUN_INIT}\n"${MPIEXEC} ${PROCLIST[i]} ${BINDIR}/${INITNAME} ${CONFLIST[i]} "|| exit 1"`
      FILES_LLIO=`echo -e ${FILES_LLIO} ${CONFLIST[i]}`
   done
fi

if [ ! ${RUNCONF} = "NONE" ]; then
   CONFLIST=(`echo ${RUNCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   FILES_LLIO=`echo -e ${FILES_LLIO} ${BINDIR}/${BINNAME}`
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      RUN_MAIN=`echo -e "${RUN_MAIN}\n"${MPIEXEC} ${PROCLIST[i]} ${BINDIR}/${BINNAME} ${CONFLIST[i]} "|| exit 1"`
      FILES_LLIO=`echo -e ${FILES_LLIO} ${CONFLIST[i]}`
   done
fi

if [ ! ${N2GCONF} = "NONE" ]; then
   CONFLIST=(`echo ${N2GCONF} | tr -s ',' ' '`)
   ndata=${#CONFLIST[@]}
   FILES_LLIO=`echo -e ${FILES_LLIO} ${BINDIR}/${N2GNAME}`
   for n in `seq 1 ${ndata}`
   do
      let i="n - 1"
      RUN_N2G=`echo -e "${RUN_N2G}\n"${MPIEXEC} ${PROCLIST[i]} ${BINDIR}/${N2GNAME} ${CONFLIST[i]} "|| exit 1"`
      FILES_LLIO=`echo -e ${FILES_LLIO} ${CONFLIST[i]}`
   done
fi

array=( `echo ${TPROC} | tr -s 'x' ' '`)
x=${array[0]}
y=${array[1]:-1}
let xy="${x} * ${y}"

if [ "${BINNAME}" = "scale-gm" ]; then
   nc=""
else
   nc=".nc"
fi

cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# ------ For FUGAKU
#
################################################################################
#PJM -g xxxxxxx
#PJM -L freq=2200
#PJM -L eco_state=2
#PJM -L rscgrp="small"
#PJM -L node=$(((TPROC+3)/4)):torus
#PJM -L elapse=01:00:00
#PJM --mpi "max-proc-per-node=4"
#PJM -j
#PJM -s
#
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#

export PARALLEL=12
export OMP_NUM_THREADS=\${PARALLEL}
export FORT90L=-Wl,-T
export PLE_MPI_STD_EMPTYFILE=off
export OMP_WAIT_POLICY=active
export FLIB_BARRIER=HARD

llio_transfer /home/apps/oss/scale/llio.list
. /home/apps/oss/scale/setup_ldlpath.sh


#. /vol0004/apps/oss/spack/share/spack/setup-env.sh
#spack load --first netcdf-c%fj
#spack load --first netcdf-fortran%fj
#spack load --first parallel-netcdf%fj

#export LD_LIBRARY_PATH=\`/home/system/tool/sort_libp\`
#/home/system/tool/sort_libp -s

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

# stage-in

llio_transfer ${FILES_LLIO}

# /home/system/tool/dir_transfer data_dir

#run

${RUN_PP}
${RUN_INIT}
${RUN_MAIN}
${RUN_N2G}


################################################################################
EOF2
