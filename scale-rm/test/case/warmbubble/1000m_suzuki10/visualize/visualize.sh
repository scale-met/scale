#! /bin/bash -x

### Energy & Mass balance ###

echo "+visualize by gnuplot"
rm -f energy.dat mass.dat mass_q.dat

while read -a line
do
   if [ ${line[0]} == "STEP=" ]; then
      echo ${line[1]} ${line[8]} ${line[9]} ${line[10]} ${line[11]} >> energy.dat
      echo ${line[1]} ${line[4]} ${line[5]} ${line[6]}  ${line[7]}  >> mass.dat
      echo ${line[1]} ${line[3]} ${line[5]} >> mass_q.dat
   fi
done < monitor.pe000000

gnuplot < ./visualize/energy.plt || exit
gnuplot < ./visualize/mass.plt   || exit
gnuplot < ./visualize/mass_q.plt || exit
rm -f energy.dat mass.dat mass_q.dat



### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

var__set=(PT W QHYD)
rangeset=(eddy auto auto)
time_set=(00000 00900 01800 02700 03600 05400)

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
   gpview history.pe\*.nc@${var},y=6000:26000,x=6000:26000,z=0:15000 --nocont --mean x,y ${eddy} --exch --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf
   gpview history.pe\*.nc@${var},y=16000,z=0:15000 --nocont --mean z ${eddy} --exch --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf column_${var}.png
   rm -f dcl.pdf

   # snapshot
   for sec in ${time_set[@]}
   do
       gpview history.pe\*.nc@${var},y=16000,z=0:15000,time=${sec} --nocont ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf ${var}${sec}sec.png
       rm -f dcl.pdf
   done

   let i="${i} + 1"
done

var__set=(RH)
rangeset=(auto auto auto auto auto auto)
time_set=

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
   gpview history.pe\*.nc@${var},y=6000:26000,x=6000:26000,z=0:15000 --nocont --mean x,y ${eddy} --exch --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf
   gpview history.pe\*.nc@${var},y=16000,z=0:15000 --nocont --mean z ${eddy} --exch --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf column_${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done
