#!/bin/sh
#-------------------------------------------------------
#   Sample Script to link parent data with grads format
#-------------------------------------------------------
rm -f ./*_?????.grd
#
dt=21600
#
# set start date (UTC)
#
start_year='2010'
start_month='05'
start_day='01'
start_hour='00'
start_min='00'
start_sec='00'
#
# set end date (UTC)
#
end_year='2010'
end_month='05'
end_day='03'
end_hour='18'
end_min='00'
end_sec='00'
#
#-----------------------------------------
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
#-----
  for BND in FNLatm FNLland FNLsfc ; do

    dir="${SCALE_DB}/FNL_output/${yyyy}${mm}"
    gpre=${BND}
    opre=${BND}

    #echo    ${dir}/${gpre}_${yyyy}${mm}${dd}${hh}.grd  ./${opre}_${fmtd_fn}.grd
    ln -svf ${dir}/${gpre}_${yyyy}${mm}${dd}${hh}.grd  ./${opre}_${fmtd_fn}.grd

  done # BND
#-----

  fn=`expr ${fn} \+ 1`
  unix_sec=`expr ${unix_sec} \+ ${dt}`
done

