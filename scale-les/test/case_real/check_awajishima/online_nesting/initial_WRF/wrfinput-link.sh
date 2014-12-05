#!/bin/sh
#-----------------------------------------
#   Test Script for XXX  (Ver.0.1)
#   2009/10/26 --- Ryuji Yoshida.
#   2014/07/08 --- Tsuyoshi Yamaura
#-----------------------------------------
dir=${SCALE_DB}'/WRF_output/kobe_case2/500m'
ftype='out'
domain='2'
#
dt=600
#
# set start date (UTC)
#
start_year='2011'
start_month='09'
start_day='20'
start_hour='03'
start_min='00'
start_sec='00'
#
# set end date (UTC)
#
end_year='2011'
end_month='09'
end_day='20'
end_hour='03'
end_min='30'
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
  year=`date -u -d "@${unix_sec}" +%Y`
  month=`date -u -d "@${unix_sec}" +%m`
  day=`date -u -d "@${unix_sec}" +%d`
  hour=`date -u -d "@${unix_sec}" +%H`
  min=`date -u -d "@${unix_sec}" +%M`
  sec=`date -u -d "@${unix_sec}" +%S`

  fmtd_fn=`printf "%05d" $fn`
  ln -svf ${dir}/wrf${ftype}_d0${domain}_${year}-${month}-${day}_${hour}\:${min}\:${sec} \
               ./wrf${ftype}_d0${domain}_${year}-${month}_${fmtd_fn}

  fn=`expr ${fn} \+ 1`
  unix_sec=`expr ${unix_sec} \+ ${dt}`
done
