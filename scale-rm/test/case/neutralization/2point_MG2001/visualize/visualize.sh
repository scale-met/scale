#! /bin/bash -x

### Energy & Mass balance ###

### Visalization ###
echo "+visualize by gpview"
rm -f dcl_*.png

var__set=(Qneut LTpath FlashPoint Epot CRGD_TOT)

i=0
for var in ${var__set[@]}
do

   # time series
   gpview history.pe\*.nc@${var},y=0.45,time=2 --nocont --wsn 2 -sw:ifl=1 || exit
   mv dcl_0001.png ${var}.png

   let i="${i} + 1"
done
gpvect --slice y=0.45,time=2 --scalar history.pe00000\*.nc@Eabs history.pe00000\*.nc@Ex history.pe00000\*.nc@Ez --xintv 2 --yintv 2 --factor 1 --unit_vect --nocont --range 0:150 --wsn 2 -sw:ifl=1 || exit
mv dcl_0001.png Evect.png
