#! /bin/bash -x

for x in 144 96 64 32 16 8 4 2 1
do
for f in 4 2
do

   let y="${x} * ${f}"
   let xy="${x} * ${y}"

   if [ ${xy} -gt 36864 ]; then
      rscgrp="huge"
   elif [ ${xy} -gt 384 ]; then
      rscgrp="large"
   else
      rscgrp="small"
   fi

   echo ${x} x ${y} = ${xy}
   cx=`printf %03d ${x}`
   cy=`printf %03d ${y}`

   expdir=005m_${x}x${y}
   mkdir -p ${expdir}

   sed -e "s/PRC_NUM_X=2/PRC_NUM_X=${x}/g" ./005m/init.conf | \
   sed -e "s/PRC_NUM_Y=3/PRC_NUM_Y=${y}/g"                  > ./${expdir}/init.conf

   sed -e "s/PRC_NUM_X=2/PRC_NUM_X=${x}/g" ./005m/run.conf  | \
   sed -e "s/PRC_NUM_Y=3/PRC_NUM_Y=${y}/g"                  > ./${expdir}/run.conf

   sed -e "s/TPROC     = 6/TPROC     = ${x}x${y}/g" ./005m/Makefile  | \
   sed -e "s/short/${rscgrp}/g"                                      > ./${expdir}/Makefile

   cp ./005m/inc_index_005m1256x32x32.f90 ./${expdir}/inc_index_005m1256x32x32.f90

   #echo "\t\$(MAKE) run -C ./005m_${x}x${y}" >> Makedef

done
done