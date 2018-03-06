#! /bin/bash -x

### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

var__set=(CN CCN)
rangeset=(auto auto)
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
   gpview history.pe\*.nc@${var},x=0,y=0,z=0 --nocont --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done
