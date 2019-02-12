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
done < monitor.peall

gnuplot < ./visualize/energy.plt || exit
gnuplot < ./visualize/mass.plt   || exit
rm -f energy.dat mass.dat



### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

var__set=(PT U W)
rangeset=(297:300 -15:15 -20:10)
time_set=(00000 00216 00396 00576 00720)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[${i}]} == "auto" ]; then
      eddy=""
      range=""
   elif [ ${rangeset[${i}]} == "eddy" ]; then
      eddy="--eddy time"
      range="--range=-3:3 --eddy x"
   else
      eddy=""
      range="--range="${rangeset[${i}]}
   fi

   # time series
   gpview history.pe\*.nc@${var},y=12000,x=9000,z=0:10000 --nocont ${eddy} --exch --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf
   gpview history.pe\*.nc@${var},y=12000,z=0 --nocont --exch --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf sfc_${var}.png
   rm -f dcl.pdf

   # snapshot
   for sec in ${time_set[@]}
   do
       gpview history.pe\*.nc@${var},y=12000,z=0:10000,time=${sec} --nocont ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf ${var}${sec}sec.png
       rm -f dcl.pdf
   done

   let i="${i} + 1"
done
