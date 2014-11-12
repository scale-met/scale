#!/bin/sh
#-----------------------------------------
#   Test Script for XXX  (Ver.0.1)
#   2009/10/26 --- Ryuji Yoshida.
#   2014/07/08 --- Tsuyoshi Yamaura
#-----------------------------------------
#
# set data dir
#
#-----------------------------------------
#idir="/data3/kenshi/jw06/g7r1_z80ST_ptb_notopo"
idir="/work/ryoshida/nicam/for_scale/jablonowski/gl07rl01pe40/g7r1_z80ST_ptb_notopo"
odir="data_source/jw06"
index='.peall'
#-----------------------------------------
#
if [ ! -d ${odir} ] ;then
 mkdir -p ${odir}
fi

ln -svf ${idir}/prs.nc    ${odir}/prs${index}.nc
ln -svf ${idir}/t.nc      ${odir}/t${index}.nc
ln -svf ${idir}/u.nc      ${odir}/u${index}.nc
ln -svf ${idir}/v.nc      ${odir}/v${index}.nc
ln -svf ${idir}/w.nc      ${odir}/w${index}.nc
ln -svf ${idir}/ps.nc     ${odir}/ps${index}.nc
  
