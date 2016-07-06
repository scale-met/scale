#!/bin/sh
for time in  600 1200  1800  2400 3000 3600
do
#gpvect --scalar --slice y=500,time=$time --noflow_vect --nocont --aspect=1.0 --xscale=20 --yscale=20 --xintv=4 --yintv=4 --unit_vect history.pe00000*@QHYD history.pe00000*@U history.pe00000*@W
gpvect --scalar --slice y=10000,time=$time --noflow_vect --nocont --aspect=1.0 --xscale=20 --yscale=20 --xintv=4 --yintv=4 --unit_vect history.pe00000*@QHYD history.pe00000*@U history.pe00000*@W
done
