#! /bin/bash -x

### Energy & Mass balance ###

echo "+visualize by gnuplot"
rm -f energy.dat mass.dat mass_q.dat

while read -a line
do
   if [ ${line[0]} == "STEP=" ]; then
      echo ${line[1]} ${line[13]} ${line[14]} ${line[15]} ${line[16]} >> energy.dat
      echo ${line[1]} ${line[9]}  ${line[10]} ${line[11]} ${line[12]} >> mass.dat
      echo ${line[1]} ${line[3]}  ${line[4]}  ${line[5]}  ${line[6]} ${line[7]} ${line[8]} ${line[10]} >> mass_q.dat
   fi
done < monitor.peall

gnuplot < ./visualize/energy.plt || exit
gnuplot < ./visualize/mass.plt   || exit
gnuplot < ./visualize/mass_q.plt || exit
rm -f energy.dat mass.dat mass_q.dat



### Visalization ###
echo "+visualize by gpview"
rm -f dcl_*.png

var__set=(QCRG_C QCRG_R QCRG_I QCRG_S QCRG_G CRGD_TOT QSPLT_G QSPLT_S QSPLT_I Qneut LTpath FlashPoint Ez Epot Eabs)

i=0
for var in ${var__set[@]}
do

   # time series
   gpview history.pe\*.nc@${var},y=6000:26000,x=6000:26000,z=0:15000 --nocont --mean x,y --exch --wsn 2 -sw:ifl=1 || exit
   mv dcl_0001.png ${var}.png

   let i="${i} + 1"
done
