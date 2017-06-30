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

var__set=(U V W PT QV)
rangeset=(auto auto auto auto auto)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time series
   gpview history.pe\*.nc@${var},z=0 --nocont --mean y --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf hov_${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done

var__set=(SFC_TEMP LHFLX SHFLX GHFLX SFLX_LW_up SFLX_LW_dn SFLX_SW_up SFLX_SW_dn U10 T2 Q2)
rangeset=(auto auto auto auto auto auto auto auto)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time series
   gpview history.pe\*.nc@${var} --nocont --mean y --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf hov_${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done

var__set=(U)
rangeset=(auto)
time_set=(1.315e+7)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # snapshot
   for sec in ${time_set[@]}
   do
       gpview history.pe\*.nc@${var},z=0:2000,time=${sec} --nocont --mean y ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf Snap_UW.${sec}sec.png
       rm -f dcl.pdf
   done

   let i="${i} + 1"
done

var__set=(U)
rangeset=(auto)
time_set=(1.316e+7 1.317e+7 1.318e+7 1.319e+7 1.320e+7)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # snapshot
   for sec in ${time_set[@]}
   do
       gpvect --scalar history.pe\*.nc@${var},z=0:2000,time=${sec} \
                       history.pe\*.nc@U,z=0:2000,time=${sec}      \
                       history.pe\*.nc@W,z=0:2000,time=${sec}      \
              --nocont --unit_vect --xintv 6 --mean y ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf Snap_UW.${sec}sec.png
       rm -f dcl.pdf
   done

   let i="${i} + 1"
done
