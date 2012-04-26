#! /bin/bash -x

tracer=ndw6

# System define
SCALEDIR=/home1/user0117/scale3

cd ${SCALEDIR}/src/project/DYCOMS2RF01

echo "##### SCALE3 auto test suite #####"
echo "SYTEM:"     ${SCALE_SYS}
echo "DIRECTORY:" ${SCALEDIR}

### 30m small job
index=030m88x8x8

echo "+++++" ${index}_${tracer}_${postfix}
echo "compile start at " `date`
make clean
make depends
make index=${index} tracer=${tracer} postfix=${postfix}

if [ ! -f ${SCALEDIR}/bin/${SCALE_SYS}/scale3_${index}_${tracer}_${postfix} ]; then
   echo "compile does not complete!" ${index}_${tracer}_${postfix}
   exit
fi

### 5m small job
index=005m438x16x16

echo "+++++" ${index}_${tracer}_${postfix}
echo "compile start at " `date`
make clean
make depends
make index=${index} tracer=${tracer} postfix=${postfix}

if [ ! -f ${SCALEDIR}/bin/${SCALE_SYS}/scale3_${index}_${tracer}_${postfix} ]; then
   echo "compile does not complete!" ${index}_${tracer}_${postfix}
   exit
fi

### 5m large job
index=005m438x8x8

echo "+++++" ${index}_${tracer}_${postfix}
echo "compile start at " `date`
make clean
make depends
make index=${index} tracer=${tracer} postfix=${postfix}

if [ ! -f ${SCALEDIR}/bin/${SCALE_SYS}/scale3_${index}_${tracer}_${postfix} ]; then
   echo "compile does not complete!" ${index}_${tracer}_${postfix}
   exit
fi

cp -r /home1/user0117/scale3/bin/K/* /data1/user0117/scale3/bin/K/