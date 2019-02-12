#! /bin/bash -x

### Energy & Mass balance ###
echo "+visualize by gnuplot"
rm -f energy.dat mass.dat mass_q.dat

while read -a line
do
   if [ ${line[0]} == "STEP=" ]; then
      echo ${line[1]} ${line[10]} ${line[11]} ${line[12]} ${line[13]} >> energy.dat
      echo ${line[1]} ${line[6]}  ${line[7]}  ${line[8]}  ${line[9]}  >> mass.dat
      echo ${line[1]} ${line[3]}  ${line[4]}  ${line[5]}  ${line[7]}  >> mass_q.dat
   fi
done < monitor.peall

gnuplot < ./visualize/energy.plt || exit
gnuplot < ./visualize/mass.plt   || exit
gnuplot < ./visualize/mass_q.plt || exit
rm -f energy.dat mass.dat mass_q.dat



### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

var__set=(U W QR)
rangeset=(auto auto 5e-6:5e-4)
time_set=(00000 00030 00060 00090)

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

   # snapshot
   for sec in ${time_set[@]}
   do
       gpview history.pe\*.nc@${var},y=500,time=${sec} --nocont ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf ${var}${sec}sec.png
       rm -f dcl.pdf
   done

   let i="${i} + 1"
done
