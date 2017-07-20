#! /bin/bash -x

### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

# initial conditions
var__set=(DENS MOMZ MOMX MOMY RHOT QV  )
rangeset=(auto auto auto auto auto auto)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   if [ ${var} == "MOMZ" ]; then
      level="zh"
   else
      level="z"
   fi

   # snapshot
   gpview init_*.pe\*.nc@${var},${level}=0 --nocont --aspect 1 ${range} --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf init.${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done

# boundary conditions
var__set=(DENS VELZ VELX VELY POTT QV  )
rangeset=(auto auto auto auto auto auto)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # average
   gpview boundary.pe\*.nc@${var},z=0 --nocont --aspect 1 --mean time ${range} --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf bnd.${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done
