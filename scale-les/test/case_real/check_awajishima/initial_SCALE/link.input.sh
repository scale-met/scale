#!/bin/sh
#-----------------------------------------
#   Test Script for XXX  (Ver.0.1)
#   2009/10/26 --- Ryuji Yoshida.
#   2014/07/08 --- Tsuyoshi Yamaura
#-----------------------------------------
dir="${SCALE_DB}/SCALE_output/japan10km_nicam"
num=9
#
#-----------------------------------------
#
fn=0

while [ ${fn} -lt ${num} ];
do
  ln -svf ${dir}/history.pe`printf "%06d" $fn`.nc .
  fn=`expr ${fn} \+ 1`
done

ln -svf ${dir}/latlon_domain_catalogue.txt .
