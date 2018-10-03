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
done < monitor.pe000000

gnuplot < ./visualize/energy.plt || exit
gnuplot < ./visualize/mass.plt   || exit
rm -f energy.dat mass.dat



### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

gpview history.pe\*.nc@W,x=42500:60000,y=50,z=0:10000,time=6000 --nozero --noshade --cint=0.2 --range=-4:4 --aspect=1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf W.png || exit

rm -f dcl.pdf
