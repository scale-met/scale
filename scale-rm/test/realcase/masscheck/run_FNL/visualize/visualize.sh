#! /bin/bash -x

### Energy & Mass balance ###

echo "+visualize by gnuplot"
rm -f energy.dat energy_flx.dat mass.dat mass_q.dat

for domain in d01 d02
do

while read -a line
do
   if [ ${line[0]} == "STEP=" ]; then
      echo ${line[1]} ${line[13]} ${line[14]} ${line[15]} ${line[16]} ${line[17]} >> energy.dat
      echo ${line[1]} ${line[17]} ${line[18]} ${line[19]} ${line[20]} ${line[21]} >> energy_flx.dat
      echo ${line[1]} ${line[9]}  ${line[10]} ${line[11]} ${line[12]}             >> mass.dat
      echo ${line[1]} ${line[3]}  ${line[4]}  ${line[5]}  ${line[6]}  ${line[7]} ${line[8]} ${line[10]} >> mass_q.dat
   fi
done < monitor_${domain}.pe000000

gnuplot < ./visualize/energy.plt     || exit
gnuplot < ./visualize/energy_flx.plt || exit
gnuplot < ./visualize/mass.plt       || exit
gnuplot < ./visualize/mass_q.plt     || exit
rm -f energy.dat energy_flx.dat mass.dat mass_q.dat

mv Energy.png     Energy_${domain}.png
mv Energy_FLX.png Energy_FLX_${domain}.png
mv Mass.png       Mass_${domain}.png
mv Mass_q.png     Mass_q_${domain}.png

done



### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

for domain in d01 d02
do

var__set=(PRES T    U    V    W    QV   QHYD)
rangeset=(auto auto auto auto auto auto auto)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # average
   gpview history_${domain}.pe\*.nc@${var},z=0 --nocont --aspect 1 --mean time --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf hist.${var}_${domain}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done

var__set=(U10  V10  T2   Q2   SFC_PRES SFC_TEMP PREC LHFLX SHFLX GHFLX SFLX_LW_dn SFLX_SW_dn)
rangeset=(auto auto auto auto auto     auto     auto auto  auto  auto  auto       auto      )

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # average
   gpview history_${domain}.pe\*.nc@${var} --nocont --aspect 1 --mean time --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf hist.${var}_${domain}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done

done
