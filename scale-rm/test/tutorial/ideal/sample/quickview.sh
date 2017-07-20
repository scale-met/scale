#!/bin/sh

echo "+visualize by gpview"
rm -f dcl.pdf

time_set=(0600 1200 1800 2400 3000 3600)

# snapshot
for sec in ${time_set[@]}
do
    gpvect --scalar --slice z=0:16000,y=0,time=${sec} --nocont --noflow_vect --unit_vect --aspect=1.0 --xscale=20 --yscale=20 --xintv=2 --yintv=2 history.pe\*.nc@QHYD history.pe\*.nc@U history.pe\*.nc@W --wsn 2 || exit
    convert -density 150 -rotate 90 +antialias dcl.pdf QHYD${sec}sec.png
    rm -f dcl.pdf
done
