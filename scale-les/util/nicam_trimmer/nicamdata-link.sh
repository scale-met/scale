#!/bin/sh
#-----------------------------------------
#   Test Script for XXX  (Ver.0.1)
#   2009/10/26 --- Ryuji Yoshida.
#   2014/07/08 --- Tsuyoshi Yamaura
#-----------------------------------------
#
# set start date (UTC)
#
start_year='2007'
start_month='10'
start_day='01'
start_hour='00'
start_min='00'
start_sec='00'
#
# set end date (UTC)
#
end_year='2007'
end_month='10'
end_day='03'
end_hour='00'
end_min='00'
end_sec='00'
#
dt=86400
#-----------------------------------------

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

#-----------------------------------------
idir="/data2/amip/${year}${month}${day}/02560x01280.zorg.torg/"
#odir="data_source/${year}${month}${day}"
odir="data_source"
index='.peall'
#-----------------------------------------
#
if [ ! -d ${odir} ] ;then
 mkdir -p ${odir}
fi

ln -svf ${idir}/la_tg.nc       ${odir}/la_tg${index}.nc
ln -svf ${idir}/la_wg.nc       ${odir}/la_wg${index}.nc
ln -svf ${idir}/ms_pres.nc     ${odir}/ms_pres${index}.nc
ln -svf ${idir}/ms_qv.nc       ${odir}/ms_qv${index}.nc
ln -svf ${idir}/ms_rh.nc       ${odir}/ms_rh${index}.nc
ln -svf ${idir}/ms_tem.nc      ${odir}/ms_tem${index}.nc
ln -svf ${idir}/ms_u.nc        ${odir}/ms_u${index}.nc
ln -svf ${idir}/ms_v.nc        ${odir}/ms_v${index}.nc
ln -svf ${idir}/oa_ice.nc      ${odir}/oa_ice${index}.nc
ln -svf ${idir}/oa_sst.nc      ${odir}/oa_sst${index}.nc
ln -svf ${idir}/sa_tppn.nc     ${odir}/sa_tppn${index}.nc
ln -svf ${idir}/ss_q2m.nc      ${odir}/ss_q2m${index}.nc
ln -svf ${idir}/ss_slp.nc      ${odir}/ss_slp${index}.nc
ln -svf ${idir}/ss_t2m.nc      ${odir}/ss_t2m${index}.nc
ln -svf ${idir}/ss_tem_sfc.nc  ${odir}/ss_tem_sfc${index}.nc
ln -svf ${idir}/ss_u10m.nc     ${odir}/ss_u10m${index}.nc
ln -svf ${idir}/ss_v10m.nc     ${odir}/ss_v10m${index}.nc

  fn=`expr ${fn} \+ 1`
  unix_sec=`expr ${unix_sec} \+ ${dt}`
done

