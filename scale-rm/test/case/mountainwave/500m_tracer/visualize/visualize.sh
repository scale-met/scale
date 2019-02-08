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

# snapshot
gpview history.pe\*.nc@W,y=500,z=0:10000,time=18000 --nozero --noshade --cint 0.1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf W.png || exit
gpview history.pe\*.nc@NC,y=500,z=0:10000,time=16200 --nocont --range=0.005:0.5 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf NC.png || exit

rm -f dcl.pdf
