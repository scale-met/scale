#! /bin/bash -x

### Energy & Mass balance ###
echo "+visualize by gnuplot"
rm -f energy.dat energy_flx.dat mass.dat energy_ocean.dat mass_ocean.dat

while read -a line
do
   if [ ${line[0]} == "STEP=" ]; then
      echo ${line[1]} ${line[5]}  ${line[6]}  ${line[7]}  ${line[8]}  ${line[9]}  >> energy.dat
      echo ${line[1]} ${line[9]}  ${line[10]} ${line[11]} ${line[12]} ${line[13]} ${line[14]} ${line[15]} >> energy_flx.dat
      echo ${line[1]} ${line[3]}  ${line[4]}                                      >> mass.dat
      echo ${line[1]} ${line[19]} ${line[20]} ${line[21]} ${line[22]} ${line[23]} ${line[24]} ${line[25]} ${line[26]} ${line[27]} >> energy_ocean.dat
      echo ${line[1]} ${line[16]} ${line[17]} ${line[18]}                         >> mass_ocean.dat
   fi
done < monitor.peall

gnuplot < ./visualize/energy.plt       || exit
gnuplot < ./visualize/energy_flx.plt   || exit
gnuplot < ./visualize/mass.plt         || exit
gnuplot < ./visualize/energy_ocean.plt || exit
gnuplot < ./visualize/mass_ocean.plt   || exit
rm -f energy.dat energy_flx.dat mass.dat energy_ocean.dat mass_ocean.dat



### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

var__set=(SFC_TEMP LHFLX SHFLX GHFLX Uabs10 T2 Q2 SFLX_LW_up SFLX_LW_dn SFLX_SW_up SFLX_SW_dn)
rangeset=(auto auto auto auto auto auto auto auto auto auto auto)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time series
   gpview history.pe\*.nc@${var},x=0,y=0 --aspect=2 --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done

var__set=(OCEAN_TEMP)
rangeset=(auto)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time series
   gpview history.pe\*.nc@${var},x=0,y=0,oz=0 --aspect=2 --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done

var__set=(U V W T)
rangeset=(auto auto auto auto)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time series
   gpview history.pe\*.nc@${var},x=0,y=0,z=0 --aspect=2 --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done

var__set=(OCEAN_ICE_TEMP OCEAN_ICE_MASS OCEAN_ICE_FRAC)
rangeset=(auto auto auto)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time series
   gpview history.pe\*.nc@${var},x=0,y=0 --aspect=2 --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done
