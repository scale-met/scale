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

var__set=(PT W TKE_SMG TKE_RS PT_W_PRIM)
rangeset=(auto auto auto auto auto)
time_set=(00000 01800 03600 05400)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[${i}]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[${i}]}
   fi

   # time series
   gpview history.pe0000\*.nc@${var} --nocont --mean x,y --exch --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf
   gpview history.pe0000\*.nc@${var} --nocont --mean y,z --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf hov_${var}.png
   rm -f dcl.pdf

   # snapshot
   for sec in ${time_set[@]}
   do
       gpview history.pe0000\*.nc@${var},y=0,time=${sec} --nocont ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf ${var}${sec}sec.png
       rm -f dcl.pdf
   done

   let i="${i} + 1"
done
