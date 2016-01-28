#!/bin/sh
#-----------------------------------------
#   Sample Script for GSM file
#-----------------------------------------

rm -f ./*_?????.grd

for BND in GSM FNL MGDSST ; do
#
dt=21600
#
# set start date (UTC)
#
start_year='2010'
start_month='06'
start_day='01'
start_hour='00'
start_min='00'
start_sec='00'
#
# set end date (UTC)
#
end_year='2010'
end_month='06'
end_day='01'
end_hour='06'
end_min='00'
end_sec='00'
#
#-----------------------------------------
case ${BND} in
GSM)
 dir=${SCALE_DB}'/GrADS_output/GSM_output'
 gpre=GANALjp
 opre=GSMjp
 ;;
FNL)
 dir=${SCALE_DB}'/GrADS_output/FNL4GSM_output'
 gpre=FNL4GANALjp
 opre=FNL
 ;;
MGDSST)
 dir=${SCALE_DB}'/GrADS_output/MGDSST_output'
 gpre=mgdsst
 opre=MGDSST
 ;;
esac
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
MGDSST)
   #echo    ${dir}/${yyyy}${mm}/${gpre}_${yyyy}${mm}${dd}${hh}.grd  ./${opre}_${fmtd_fn}.grd
   ln -svf ${dir}/${yyyy}/${gpre}_${yyyy}${mm}${dd}.grd  ./${opre}_${fmtd_fn}.grd
;;
*)
   #echo    ${dir}/${yyyy}${mm}/${gpre}_${yyyy}${mm}${dd}${hh}.grd  ./${opre}_${fmtd_fn}.grd
   ln -svf ${dir}/${yyyy}${mm}/${gpre}_${yyyy}${mm}${dd}${hh}.grd  ./${opre}_${fmtd_fn}.grd
;;
esac

  fn=`expr ${fn} \+ 1`
  unix_sec=`expr ${unix_sec} \+ ${dt}`
done

done # BND
