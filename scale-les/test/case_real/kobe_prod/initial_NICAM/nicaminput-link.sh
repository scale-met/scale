#!/bin/sh
#-----------------------------------------
#   Test Script for XXX  (Ver.0.1)
#   2009/10/26 --- Ryuji Yoshida.
#   2014/07/08 --- Tsuyoshi Yamaura
#-----------------------------------------
dir='/work3/scale/external_data/nicam/amip/20071001'
index='.peall'
#-----------------------------------------
#

ln -svf ${dir}/la_tg.nc       ./la_tg${index}.nc
ln -svf ${dir}/la_wg.nc       ./la_wg${index}.nc
ln -svf ${dir}/ms_pres.nc     ./ms_pres${index}.nc
ln -svf ${dir}/ms_qv.nc       ./ms_qv${index}.nc
ln -svf ${dir}/ms_rh.nc       ./ms_rh${index}.nc
ln -svf ${dir}/ms_tem.nc      ./ms_tem${index}.nc
ln -svf ${dir}/ms_u.nc        ./ms_u${index}.nc
ln -svf ${dir}/ms_v.nc        ./ms_v${index}.nc
ln -svf ${dir}/oa_ice.nc      ./oa_ice${index}.nc
ln -svf ${dir}/oa_sst.nc      ./oa_sst${index}.nc
ln -svf ${dir}/sa_tppn.nc     ./sa_tppn${index}.nc
ln -svf ${dir}/ss_q2m.nc      ./ss_q2m${index}.nc
ln -svf ${dir}/ss_slp.nc      ./ss_slp${index}.nc
ln -svf ${dir}/ss_t2m.nc      ./ss_t2m${index}.nc
ln -svf ${dir}/ss_tem_sfc.nc  ./ss_tem_sfc${index}.nc
ln -svf ${dir}/ss_u10m.nc     ./ss_u10m${index}.nc
ln -svf ${dir}/ss_v10m.nc     ./ss_v10m${index}.nc

