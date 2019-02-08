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

var__set=(W WA WZ WP MOMZ MOMZ_Z MOMZ_P)
rangeset=
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
   gpview history.pe\*.nc@${var},x=30000,y=20000 --nocont --exch --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf

   # snapshot
   for sec in ${time_set[@]}
   do
       gpview history.pe\*.nc@${var},time=${sec} --nocont --mean y ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf ${var}${sec}sec.png
       rm -f dcl.pdf
   done

   let i="${i} + 1"
done

var__set=(MOMX MOMX_Z MOMX_P)
rangeset=
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
   gpview history.pe\*.nc@${var},xh=30000,y=20000 --nocont --exch --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf

   # snapshot
   for sec in ${time_set[@]}
   do
       gpview history.pe\*.nc@${var},time=${sec} --nocont --mean y ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf ${var}${sec}sec.png
       rm -f dcl.pdf
   done

   let i="${i} + 1"
done

var__set=(MOMY MOMY_Z MOMY_P)
rangeset=
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
   gpview history.pe\*.nc@${var},x=30000,yh=20000 --nocont --exch --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf

   # snapshot
   for sec in ${time_set[@]}
   do
       gpview history.pe\*.nc@${var},time=${sec} --nocont --mean y ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf ${var}${sec}sec.png
       rm -f dcl.pdf
   done

   let i="${i} + 1"
done
