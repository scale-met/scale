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
MPIEXEC="mpirun -np ${TPROC}"

# Generate run.sh

cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# ------ FOR Linux64 & intel C&fortran & intel mpi -----
#
################################################################################
export FORT_FMT_RECL=400

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
   echo "${MPIEXEC} ${BINDIR}/${BINNAME}  ${RUNCONF}  || exit 1" >> ./run.sh
fi
