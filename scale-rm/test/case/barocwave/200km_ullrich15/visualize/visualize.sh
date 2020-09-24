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

var__set=(PRES T)
rangeset=(910e2:958e2 264:304)
intset=(2e2 2)
time_set=(00000 518400 691200 864000 1036800)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   if [ ${intset[$i]} == "auto" ]; then
      intset=""
   else
      intset="--int="${intset[$i]}
   fi

   # time series
   #gpview history.pe\*.nc@${var},x=0,y=0,z=5000,time=400:500 --wsn 2 || exit
   #convert -density 150 -rotate 90 +antialias dcl.pdf slice_${var}.png
   #rm -f dcl.pdf

   # snapshot
   for sec in ${time_set[@]}
   do
       gpview history.pe\*.nc@${var},time=${sec},z=500,x=0:4e7 --nocont ${range} ${intset} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf ${var}${sec}sec.png
       rm -f dcl.pdf
   done

   let i="${i} + 1"
done
