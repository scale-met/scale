#! /bin/bash -x

### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

var__set=(DENS RHOT QV)
rangeset=(auto auto auto)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # restart 0
   gpview restartA_\*.pe\*.nc@${var},z=0 --nocont --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf RestartA_${var}.png
   rm -f dcl.pdf
   # restart B
   gpview restartB2_\*.pe\*.nc@${var},z=0 --nocont --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf RestartB2_${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done
