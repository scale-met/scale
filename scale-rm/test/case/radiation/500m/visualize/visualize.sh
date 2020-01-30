#! /bin/bash -x

### Energy & Mass balance ###
echo "+visualize by gnuplot"
rm -f energy.dat energy_sfc.dat energy_tom.dat

while read -a line
do
   if [ ${line[0]} == "STEP=" ]; then
      echo ${line[1]} ${line[5]} ${line[6]}  ${line[7]}                          >> energy.dat
      echo ${line[1]} ${line[6]} ${line[8]}  ${line[9]}  ${line[10]} ${line[11]} >> energy_sfc.dat
      echo ${line[1]} ${line[7]} ${line[12]} ${line[13]} ${line[14]} ${line[15]} >> energy_tom.dat
   fi
done < monitor.peall

gnuplot < ./visualize/energy.plt || exit
gnuplot < ./visualize/energy_sfc.plt || exit
gnuplot < ./visualize/energy_tom.plt || exit
rm -f energy.dat energy_sfc.dat energy_tom.dat



### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

var__set=(TEMP_t_rd TEMP_t_rd_LW TEMP_t_rd_SW RADFLUX_LW RADFLUX_SW)
rangeset=(auto auto auto auto auto)
time_set=

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time series
   gpview history_\*.pe\*.nc@${var},x=0,y=0 --nocont --exch --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf

   # snapshot
   for sec in ${time_set[@]}
   do
       gpview history_\*.pe\*.nc@${var},time=${sec} --nocont --mean y ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf ${var}${sec}sec.png
       rm -f dcl.pdf
   done

   let i="${i} + 1"
done

rm -f dcl.pdf

var__set=(COSZ SOLINS OLR OSR SLR SSR)
rangeset=(auto auto auto auto auto auto)
time_set=

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time series
   gpview history_\*.pe\*.nc@${var},x=0,y=0 --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf

   # snapshot
   for sec in ${time_set[@]}
   do
       gpview history_\*.pe\*.nc@${var},time=${sec} --nocont --mean y ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf ${var}${sec}sec.png
       rm -f dcl.pdf
   done

   let i="${i} + 1"
done
