#! /bin/bash -x

# Arguments
BINDIR=${1}
UTILDIR=${2}
PPNAME=${3}
INITNAME=${4}
BINNAME=${5}
N2GNAME=${6}
PPCONF=${7}
INITCONF=${8}
RUNCONF=${9}
N2GCONF=${10}
TPROC=${11}
DATDIR=${12}
DATPARAM=(`echo ${13} | tr -s ',' ' '`)
DATDISTS=(`echo ${14} | tr -s ',' ' '`)

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

if [ ! ${N2GCONF} = "NONE" ]; then
  RUN_N2G="${MPIEXEC} ${UTILDIR}/${N2GNAME} ${N2GCONF} || exit"
fi





cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# ------ For Linux64 & gnu fortran&C & openmpi -----
#
################################################################################
export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

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
${RUN_N2G}

################################################################################
EOF2
