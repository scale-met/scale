#! /bin/bash -x

for domain in d01 d02
do

### Energy & Mass balance ###

echo "+visualize by gnuplot"
rm -f energy.dat energy_flx.dat mass.dat mass_q.dat

while read -a line
do
   if [ ${line[0]} == "STEP=" ]; then
      echo ${line[1]} ${line[7]}  ${line[8]}  ${line[9]}  ${line[10]} ${line[11]} >> energy.dat
      echo ${line[1]} ${line[11]} ${line[12]} ${line[13]} ${line[14]} ${line[15]} >> energy_flx.dat
      echo ${line[1]} ${line[3]}  ${line[4]}  ${line[5]}  ${line[6]}              >> mass.dat
      echo ${line[1]} ${line[4]}  ${line[16]} ${line[17]} ${line[18]} ${line[19]} ${line[20]} ${line[21]} >> mass_q.dat
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
