#!/bin/sh
#-----------------------------------------
#   Test Script for XXX  (Ver.0.1)
#   2009/10/26 --- Ryuji Yoshida.
#   2014/07/08 --- Tsuyoshi Yamaura
#-----------------------------------------
dir=${SCALE_DB}'/NICAM_output/trimmed_data'
ftype='peall'
#
dt=86400
#
# set start date (UTC)
#
start_year='1999'
start_month='05'
start_day='01'
start_hour='06'
start_min='00'
start_sec='00'
#
# set end date (UTC)
#
end_year='1999'
end_month='05'
end_day='02'
end_hour='00'
end_min='00'
end_sec='00'
#
#-----------------------------------------
#
fn=0
rm -f ./*_?????.${ftype}.nc

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

  ln -svf ${dir}/${year}${month}${day}/la_tg.nc      ./la_tg_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/${year}${month}${day}/la_wg.nc      ./la_wg_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/${year}${month}${day}/lsmask.nc     ./lsmask_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/${year}${month}${day}/ms_pres.nc    ./ms_pres_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/${year}${month}${day}/ms_qv.nc      ./ms_qv_${fmtd_fn}.${ftype}.nc
  #ln -svf ${dir}/${year}${month}${day}/ms_rh.nc      ./ms_rh_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/${year}${month}${day}/ms_tem.nc     ./ms_tem_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/${year}${month}${day}/ms_u.nc       ./ms_u_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/${year}${month}${day}/ms_v.nc       ./ms_v_${fmtd_fn}.${ftype}.nc
  #ln -svf ${dir}/${year}${month}${day}/oa_ice.nc     ./oa_ice_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/${year}${month}${day}/oa_sst.nc     ./oa_sst_${fmtd_fn}.${ftype}.nc
  #ln -svf ${dir}/${year}${month}${day}/sa_tppn.nc    ./sa_tppn_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/${year}${month}${day}/ss_q2m.nc     ./ss_q2m_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/${year}${month}${day}/ss_slp.nc     ./ss_slp_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/${year}${month}${day}/ss_t2m.nc     ./ss_t2m_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/${year}${month}${day}/ss_tem_sfc.nc ./ss_tem_sfc_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/${year}${month}${day}/ss_u10m.nc    ./ss_u10m_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/${year}${month}${day}/ss_v10m.nc    ./ss_v10m_${fmtd_fn}.${ftype}.nc

  fn=`expr ${fn} \+ 1`
  unix_sec=`expr ${unix_sec} \+ ${dt}`
done
