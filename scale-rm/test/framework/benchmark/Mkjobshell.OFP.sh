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
DATDIR=${9}
DATPARAM=(`echo ${10} | tr -s ',' ' '`)
DATDISTS=(`echo ${11} | tr -s ',' ' '`)

# System specific
MPIEXEC="mpirun -np ${TPROC}"

if [ ! ${PPCONF} = "NONE" ]; then
  RUN_PP="${MPIEXEC} ${BINDIR}/${PPNAME} ${PPCONF} || exit"
fi

if [ ! ${INITCONF} = "NONE" ]; then
  RUN_INIT="${MPIEXEC} ${BINDIR}/${INITNAME} ${INITCONF} || exit"
fi

if [ ! ${RUNCONF} = "NONE" ]; then
  RUN_BIN="${MPIEXEC} ${BINDIR}/${BINNAME} ${RUNCONF} || exit"
fi

NNODE=`expr \( $TPROC - 1 \) / 64 + 1`
NPROC=`expr $TPROC / $NNODE`
NPIN=`expr 255 / \( $NPROC \) + 1`

cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# ------ FOR Oakforest-PACS -----
#
################################################################################
#PJM -g gc26
#PJM -L rscgrp=regular-cache
#PJM -L node=${NNODE}
#PJM --mpi proc=${TPROC}
#PJM --omp thread=1
#PJM -L elapse=00:30:00
#PJM -N SCALE
#PJM -X
#PJM -j
#PJM -s
#
module load hdf5_szip
module load hdf5
module load netcdf
module load netcdf-fortran

export FORT_FMT_RECL=400

export HFI_NO_CPUAFFINITY=1
export I_MPI_PIN_PROCESSOR_EXCLUDE_LIST=0,1,68,69,136,137,204,205
export I_MPI_HBW_POLICY=hbw_preferred,,
export I_MPI_FABRICS_LIST=tmi
unset KMP_AFFINITY
#export KMP_AFFINITY=verbose
#export I_MPI_DEBUG=5

export OMP_NUM_THREADS=1
export I_MPI_PIN_DOMAIN=${NPIN}
export I_MPI_PERHOST=${NPROC}
export KMP_HW_SUBSET=1t

EOF1

if [ ! ${DATPARAM[0]} = "" ]; then
   for f in ${DATPARAM[@]}
   do
         if [ -f ${DATDIR}/${f} ]; then
            echo "ln -svf ${DATDIR}/${f} ." >> ./run.sh
         elif [ -d ${DATDIR}/${f} ]; then
            echo "rm -f                  ./input" >> ./run.sh
            echo "ln -svf ${DATDIR}/${f} ./input" >> ./run.sh
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
         if [ -f ${f}.pe${PE}.nc ]; then
            echo "ln -svf ${f}.pe${PE}.nc ." >> ./run.sh
         else
            echo "datafile does not found! : ${f}.pe${PE}.nc"
            exit 1
         fi
      done
   done
fi

cat << EOF2 >> ./run.sh

# run
${RUN_PP}
${RUN_INIT}
${RUN_BIN}

################################################################################
EOF2
