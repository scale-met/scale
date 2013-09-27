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
gpview ${DATADIR}/history.pe\*.nc@W,x=95000:120000,y=100,z=0:15000,time=6000 --nozero --noshade --cint=0.1 --aspect=1 --wsn 2 || exit

echo "+convert to png"
convert -density 150 -rotate 90 +antialias dcl.ps W.png || exit
rm -f dcl.ps
