#! /bin/bash -x

### Energy & Mass balance ###
echo "+visualize by gnuplot"
rm -f energy.dat energy_flx.dat mass.dat

while read -a line
do
   if [ ${line[0]} == "STEP=" ]; then
      echo ${line[1]} ${line[5]} ${line[6]}  ${line[7]}  ${line[8]}  ${line[9]}  >> energy.dat
      echo ${line[1]} ${line[9]} ${line[10]} ${line[11]} ${line[12]} ${line[13]} >> energy_flx.dat
      echo ${line[1]} ${line[3]} ${line[4]}                                      >> mass.dat
   fi
done < monitor.pe000000

gnuplot < ./visualize/energy.plt     || exit
gnuplot < ./visualize/energy_flx.plt || exit
gnuplot < ./visualize/mass.plt       || exit
rm -f energy.dat energy_flx.dat mass.dat



### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

var__set=(LAND_TEMP LAND_WATER)
rangeset=(auto auto)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time series
   gpview history.pe\*.nc@${var},x=0,y=0 --aspect=2 --nocont --exch --itr=2 --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done

var__set=(SFC_TEMP LHFLX SHFLX GHFLX Uabs10 T2 Q2 SFLX_LW_up SFLX_LW_dn SFLX_SW_up SFLX_SW_dn)
rangeset=(auto auto auto auto 0.1:0.35 auto auto auto auto auto auto)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range=${rangeset[$i]}"
   fi

   # time series
   gpview history.pe\*.nc@${var},x=0,y=0 --aspect=2 --wsn 2 ${range} || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done
