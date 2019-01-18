#! /bin/bash -x

### Energy & Mass balance ###
echo "+visualize by gnuplot"
rm -f energy.dat mass.dat

while read -a line
do
   if [ ${line[0]} == "STEP=" ]; then
      echo ${line[1]} ${line[5]} ${line[6]} ${line[7]} ${line[8]} >> energy.dat
      echo ${line[1]} ${line[3]} ${line[4]}                       >> mass.dat
   fi
done < monitor.pe000000

gnuplot < ./visualize/energy.plt || exit
gnuplot < ./visualize/mass.plt   || exit
rm -f energy.dat mass.dat



### Visualization ###
echo "+visualize by gpview"
rm -f dcl.pdf

var__set=(SFC_TEMP SHFLX LHFLX)
rangeset=(auto auto auto)
time_set=

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time seriese
   gpview history.pe\*.nc@${var},x=0,y=0 ${range} --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done

var__set=(PT U V W QV TKE_MYNN Nu_BL Kh_BL Pr_MYNN Ri_MYNN L_mix_MYNN TKE_prod_MYNN TKE_diss_MYNN)
zcordset=(z  z z z z  z        zh    zh    zh      z       z          z             z            )
rangeset=(275:310 -1:7 -13:0 -5e-4:5e-4 1e-4:3e-3 0:2 0:150 0:150 0:10 -20:20 0:100 -2e-3:1e-2 -1e-2:0)
time_set=

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # hovmoller diagram
   gpview history.pe\*.nc@${var},x=0,y=0,${zcordset[$i]}=0:3400 --nocont ${range} --wsn 2 --exch || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf hov_${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done
