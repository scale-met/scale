#! /bin/bash -x

# export postfix=${1:-20120413}
tracer=ndw6

# System define
SCALEDIR=/home/c18001/scale3
LOGFILE=${SCALEDIR}/run/${SCALE_SYS}/compilelog_${postfix}.log

cd ${SCALEDIR}/src

echo "##### SCALE3 auto test suite #####" >  ${LOGFILE}
echo "SYTEM:" ${SCALE_SYS}                >> ${LOGFILE}
echo "DIRECTORY:" ${SCALEDIR}             >> ${LOGFILE}

### 500m experiment
index=42x56x56
echo "+++++" ${index}_${tracer}_${postfix}                  >> ${LOGFILE}
echo "compile start at " `date`                             >> ${LOGFILE}
make clean                                                  >> ${LOGFILE}
make -j2 index=${index} tracer=${tracer} postfix=${postfix} >> ${LOGFILE}

if [ ! -f ${SCALEDIR}/bin/${SCALE_SYS}/scale3_${index}_${tracer}_${postfix} ]; then
   echo "compile does not complete!" ${SCALEDIR}/bin/${SCALE_SYS}/scale3_${index}_${tracer}_${postfix}
   exit
fi

### 50m experiment
index=206x21x21

echo "+++++" ${index}_${tracer}_${postfix}                  >> ${LOGFILE}
echo "compile start at " `date`                             >> ${LOGFILE}
make clean                                                  >> ${LOGFILE}
make -j2 index=${index} tracer=${tracer} postfix=${postfix} >> ${LOGFILE}

if [ ! -f ${SCALEDIR}/bin/${SCALE_SYS}/scale3_${index}_${tracer}_${postfix} ]; then
   echo "compile does not complete!" ${index}_${tracer}_${postfix}
   exit
fi

### 5m experiment
index=436x14x14

echo "+++++" ${index}_${tracer}_${postfix}                  >> ${LOGFILE}
echo "compile start at " `date`                             >> ${LOGFILE}
make clean                                                  >> ${LOGFILE}
make -j2 index=${index} tracer=${tracer} postfix=${postfix} >> ${LOGFILE}

if [ ! -f ${SCALEDIR}/bin/${SCALE_SYS}/scale3_${index}_${tracer}_${postfix} ]; then
   echo "compile does not complete!" ${index}_${tracer}_${postfix}
   exit
fi

### 5m performance
index=1256x28x28

echo "+++++" ${index}_${tracer}_${postfix}                  >> ${LOGFILE}
echo "compile start at " `date`                             >> ${LOGFILE}
make clean                                                  >> ${LOGFILE}
make -j2 index=${index} tracer=${tracer} postfix=${postfix} >> ${LOGFILE}

if [ ! -f ${SCALEDIR}/bin/${SCALE_SYS}/scale3_${index}_${tracer}_${postfix} ]; then
   echo "compile does not complete!" ${index}_${tracer}_${postfix}
   exit
fi


### 50m narrow experiment
#index=72x176x8
#
#echo "+++++" ${index}_${tracer}_${postfix}                  >> ${LOGFILE}
#echo "compile start at " `date`                             >> ${LOGFILE}
#make clean                                                  >> ${LOGFILE}
#make -j2 index=${index} tracer=${tracer} postfix=${postfix} >> ${LOGFILE}
#
#if [ ! -f ${SCALEDIR}/bin/${SCALE_SYS}/scale3_${index}_${tracer}_${postfix} ]; then
#   echo "compile does not complete!" ${index}_${tracer}_${postfix}
#   exit
#fi
