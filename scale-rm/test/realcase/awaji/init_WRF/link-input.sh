#!/bin/sh
################################################################################
# 2009/10/26 --- Ryuji Yoshida.
# 2014/07/08 --- Tsuyoshi Yamaura
################################################################################

# set start date (UTC)
yr_str='2011'
mo_str='09'
dy_str='07'
hr_str='00'
mn_str='00'
sc_str='00'

# set _end date (UTC)
yr_end='2011'
mo_end='09'
dy_end='07'
hr_end='06'
mn_end='00'
sc_end='00'

usec_str=`date -u -d "${yr_str}-${mo_str}-${dy_str} ${hr_str}:${mn_str}:${sc_str}" +%s`
usec_end=`date -u -d "${yr_end}-${mo_end}-${dy_end} ${hr_end}:${mn_end}:${sc_end}" +%s`

################################################################################

indir="${SCALE_DB}/init_sample/WRF_2.5km"
ftype='out'
domain='1'
fname="wrf${ftype}_d0${domain}"

dt=600

fnum=0
usec=${usec_str}
while [ ${usec} -le ${usec_end} ];
do
   FN=`printf %05d ${fnum}`

   date_now=`date -u -d "@${usec}" +%Y-%m-%d_%H:%M:%S`

   if [ "${SCALE_SYS}" == "Kmicro" ]; then
      cp -vf  ${indir}/${fname}_${date_now} ./${fname}_${FN}
   else
      ln -svf ${indir}/${fname}_${date_now} ./${fname}_${FN}
   fi

   let fnum="${fnum}+1"
   let usec="${usec}+${dt}"
done
