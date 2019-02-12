#! /bin/bash -x

### Energy & Mass balance ###
echo "+visualize by gnuplot"
rm -f energy.dat mass.dat

while read -a line
do
   if [ ${line[0]} == "STEP=" ]; then
      echo ${line[1]} ${line[4]} ${line[5]} ${line[6]} ${line[7]} >> energy.dat
      echo ${line[1]} ${line[4]}                                  >> mass.dat
   fi
done < monitor.peall

gnuplot < ./visualize/energy.plt || exit
gnuplot < ./visualize/mass.plt   || exit
rm -f energy.dat mass.dat



### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

var__set=(U V PRES)
rangeset=(-0.05:0.05 -0.08:0.08 auto)
time_set=(00000 172800)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
       range=""
   else
       range="--range="${rangeset[$i]}
   fi
   if [ ${var} == "PRES" ]; then
       eddy="--eddy x"
   else
       eddy=""
   fi

   # time series
   gpview history.pe\*.nc@${var},y=4e6,z=0 --nocont --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   rm -f dcl.pdf

   # snapshot
   for sec in ${time_set[@]}
   do
       gpview history.pe\*.nc@${var},z=0,time=${sec} --nocont ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf ${var}_horiz_${sec}sec.png
       rm -f dcl.pdf

       gpview history.pe\*.nc@${var},y=4e6,time=${sec} --nocont ${range} ${eddy} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf ${var}_vert_${sec}sec.png
       rm -f dcl.pdf
   done

   let i="${i} + 1"
done
