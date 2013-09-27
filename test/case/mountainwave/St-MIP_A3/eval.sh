#! /bin/bash -x

DATADIR=${1:-"."}



### Energy & Mass balance ###
echo "+visualize by gnuplot"
rm -f energy.dat mass.dat

while read -a line
do
   if [ ${line[0]} == "STEP=" ]; then
      echo ${line[1]} ${line[5]} ${line[6]} ${line[7]} ${line[8]} >> energy.dat
      echo ${line[1]} ${line[3]} ${line[4]}                       >> mass.dat
   fi
done < ${DATADIR}/monitor.pe000000

gnuplot < ../St-MIP_A1/energy.plt || exit
gnuplot < ../St-MIP_A1/mass.plt   || exit
rm -f energy.dat mass.dat



### Visalization ###
echo "+visualize by gpview"
rm -f dcl.ps

echo "+visualize by gpview"
gpview ${DATADIR}/history.pe\*.nc@W,x=19000:21000,y=20,z=0:2000,time=1200 --nozero --noshade --cint=0.25 --aspect=1 --wsn 2 || exit

echo "+convert to png"
convert -density 150 -rotate 90 +antialias dcl.ps W.png || exit
rm -f dcl.ps
