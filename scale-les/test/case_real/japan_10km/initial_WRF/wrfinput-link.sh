#!/bin/sh
#-----------------------------------------
#   Test Script for XXX  (Ver.0.1)
#   2009/10/26 --- Ryuji Yoshida.
#   2014/07/08 --- Tsuyoshi Yamaura
#-----------------------------------------
dir="/data0/scale_database/WRF_output/14km_exp/FNL.10hPa"
ftype='out'
domain='1'
#
dt=21600
#
# set start date (UTC)
#
start_year='2011'
start_month='09'
start_day='19'
start_hour='18'
start_min='00'
start_sec='00'
#
# set end date (UTC)
#
end_year='2011'
end_month='09'
end_day='22'
end_hour='00'
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

