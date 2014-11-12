#!/bin/sh
#-----------------------------------------
#   Test Script for XXX  (Ver.0.1)
#   2009/10/26 --- Ryuji Yoshida.
#   2014/07/08 --- Tsuyoshi Yamaura
#-----------------------------------------
dir='/work/ryoshida/scale/test/dev1111/scale-les/util/nicam_trimmer/jw06'
ftype='peall'
#
#-----------------------------------------
#
fn=0
rm -f ./*_?????.${ftype}.nc

  fmtd_fn=`printf "%05d" $fn`

  ln -svf ${dir}/la_tg.nc      ./la_tg_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/la_wg.nc      ./la_wg_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/ms_pres.nc    ./ms_pres_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/ms_qv.nc      ./ms_qv_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/ms_tem.nc     ./ms_tem_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/ms_u.nc       ./ms_u_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/ms_v.nc       ./ms_v_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/oa_sst.nc     ./oa_sst_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/ss_q2m.nc     ./ss_q2m_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/ss_slp.nc     ./ss_slp_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/ss_t2m.nc     ./ss_t2m_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/ss_tem_sfc.nc ./ss_tem_sfc_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/ss_u10m.nc    ./ss_u10m_${fmtd_fn}.${ftype}.nc
  ln -svf ${dir}/ss_v10m.nc    ./ss_v10m_${fmtd_fn}.${ftype}.nc

