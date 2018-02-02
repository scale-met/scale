#!/bin/sh
################################################################################
# 2009/10/26 --- Ryuji Yoshida.
# 2014/07/08 --- Tsuyoshi Yamaura
################################################################################

# set start date (UTC)
yr_str='2000'
mo_str='07'
dy_str='01'
hr_str='00'
mn_str='00'
sc_str='00'

# set _end date (UTC)
yr_end='2000'
mo_end='07'
dy_end='04'
hr_end='00'
mn_end='00'
sc_end='00'

usec_str=`date -u -d "${yr_str}-${mo_str}-${dy_str} ${hr_str}:${mn_str}:${sc_str}" +%s`
usec_end=`date -u -d "${yr_end}-${mo_end}-${dy_end} ${hr_end}:${mn_end}:${sc_end}" +%s`

################################################################################

indir="${SCALE_DB}/init_sample/MRIAGCM_20km"

dt=21600

fnum=0
usec=${usec_str}
while [ ${usec} -le ${usec_end} ];
do
   FN=`printf %05d ${fnum}`

   year_dir=`date -u -d "@${usec}" +%Y`
   date_now=`date -u -d "@${usec}" +%Y%m%d%H`

   if [ "${SCALE_SYS}" == "Kmicro" ]; then
      cp -vf  ${indir}/invariant_TL959.dr                      ./MRI_${FN}.grd
      cp -vf  ${indir}/${year_dir}/atm_eta_pres_${date_now}.dr ./MRIpres_${FN}.grd
      cp -vf  ${indir}/${year_dir}/atm_eta_${date_now}.dr      ./MRIatm_${FN}.grd
      cp -vf  ${indir}/${year_dir}/sfc_${date_now}.dr          ./MRIsfc_${FN}.grd
   else
      ln -svf ${indir}/invariant_TL959.dr                      ./MRI_${FN}.grd
      ln -svf ${indir}/${year_dir}/atm_eta_pres_${date_now}.dr ./MRIpres_${FN}.grd
      ln -svf ${indir}/${year_dir}/atm_eta_${date_now}.dr      ./MRIatm_${FN}.grd
      ln -svf ${indir}/${year_dir}/sfc_${date_now}.dr          ./MRIsfc_${FN}.grd
   fi

   let fnum="${fnum}+1"
   let usec="${usec}+${dt}"
done
