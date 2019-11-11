#!/bin/sh
################################################################################
# 2009/10/26 --- Ryuji Yoshida.
# 2014/07/08 --- Tsuyoshi Yamaura
################################################################################

# set start date (UTC)
yr_str='2010'
mo_str='05'
dy_str='01'
hr_str='00'
mn_str='00'
sc_str='00'

# set _end date (UTC)
yr_end='2010'
mo_end='05'
dy_end='04'
hr_end='00'
mn_end='00'
sc_end='00'

usec_str=`date -u -d "${yr_str}-${mo_str}-${dy_str} ${hr_str}:${mn_str}:${sc_str}" +%s`
usec_end=`date -u -d "${yr_end}-${mo_end}-${dy_end} ${hr_end}:${mn_end}:${sc_end}" +%s`

################################################################################

indir="${SCALE_DB}/init_sample/FNL_1deg"

dt=21600

fnum=0
usec=${usec_str}
while [ ${usec} -le ${usec_end} ];
do
   FN=`printf %05d ${fnum}`

   year_dir=`date -u -d "@${usec}" +%Y`
   date_now=`date -u -d "@${usec}" +%Y%m%d%H`

   ln -svf ${indir}/${year_dir}/FNL_ATM_${date_now}.grd ./ATM_${FN}.grd
   ln -svf ${indir}/${year_dir}/FNL_SFC_${date_now}.grd ./SFC_${FN}.grd
   ln -svf ${indir}/${year_dir}/FNL_LND_${date_now}.grd ./LND_${FN}.grd

   let fnum="${fnum}+1"
   let usec="${usec}+${dt}"
done
