#!/bin/sh
#-----------------------------------------
#   Sample Script for MRI-AGCM file
#-----------------------------------------

rm -f ./*_?????.grd

for BND in MRI MRIpres MRIatm MRIsfc ; do
#
dt=21600
#
# set start date (UTC)
#
start_year='2000'
start_month='07'
start_day='01'
start_hour='00'
start_min='00'
start_sec='00'
#
# set end date (UTC)
#
end_year='2000'
end_month='07'
end_day='03'
end_hour='18'
end_min='00'
end_sec='00'
#
#-----------------------------------------
#
fn=0

start_unix_sec=`date -u -d "${start_year}-${start_month}-${start_day} ${start_hour}:${start_min}:${start_sec}" +%s`
end_unix_sec=`date -u -d "${end_year}-${end_month}-${end_day} ${end_hour}:${end_min}:${end_sec}" +%s`

unix_sec=${start_unix_sec}
while [ ${unix_sec} -le ${end_unix_sec} ];
do
  yyyy=`date -u -d "@${unix_sec}" +%Y`
  mm=`date -u -d "@${unix_sec}" +%m`
  dd=`date -u -d "@${unix_sec}" +%d`
  hh=`date -u -d "@${unix_sec}" +%H`
  nn=`date -u -d "@${unix_sec}" +%M`
  ss=`date -u -d "@${unix_sec}" +%S`

  fmtd_fn=`printf "%05d" $fn`

case ${BND} in
MRI)
 dir=${SCALE_DB}'/GrADS_output/MRIAGCM_output/regional'
 gfile=invariant_105E-170E_10N-65N_TL959.dr
 opre=${BND}
 ;;
MRIpres)
 dir=${SCALE_DB}"/GrADS_output/MRIAGCM_output/${yyyy}${mm}"
 gfile=atm_eta_pres_${yyyy}${mm}${dd}${hh}.dr
 opre=${BND}
 ;;
MRIatm)
 dir=${SCALE_DB}"/GrADS_output/MRIAGCM_output/${yyyy}${mm}"
 gfile=atm_eta_${yyyy}${mm}${dd}${hh}.dr
 opre=${BND}
 ;;
MRIsfc)
 dir=${SCALE_DB}"/GrADS_output/MRIAGCM_output/${yyyy}${mm}"
 gfile=sfc_${yyyy}${mm}${dd}${hh}.dr
 opre=${BND}
 ;;
esac
#

   #echo    ${dir}/${yyyy}${mm}/${gpre}_${yyyy}${mm}${dd}${hh}.grd  ./${opre}_${fmtd_fn}.grd
   ln -svf ${dir}/${gfile}   ./${opre}_${fmtd_fn}.grd


  fn=`expr ${fn} \+ 1`
  unix_sec=`expr ${unix_sec} \+ ${dt}`
done

done # BND
