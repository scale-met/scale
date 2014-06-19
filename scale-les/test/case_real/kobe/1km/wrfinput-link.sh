#!/bin/sh
#-----------------------------------------
#   Test Script for XXX  (Ver.0.1)
#   2009/10/26 --- Ryuji Yoshida.
#-----------------------------------------
dir='/data3/kenshi/140507_WRFcase-exp_for_scale/test5/WRF_recca/case1'
ftype='out'
domain='2'
yyyy='2011'
mm='09'
dd1='02'

ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:00\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00000
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:02\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00001
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:04\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00002
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:06\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00003
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:08\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00004
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:10\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00005
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:12\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00006
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:14\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00007
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:16\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00008
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:18\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00009
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:20\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00010
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:22\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00011
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:24\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00012
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:26\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00013
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:28\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00014
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:30\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00015
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:32\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00016
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:34\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00017
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:36\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00018
ln -s ${dir}/wrf${ftype}_d0${domain}_${yyyy}-${mm}-${dd1}_06\:38\:00 ./wrf${ftype}_d0${domain}_${yyyy}-${mm}_00019

#-----------------------------------------
#EOF
