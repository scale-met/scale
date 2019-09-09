#/bin/bash

gpview --exch --range 0:0.018 --overplot 2 ../history.pe00000\*.nc@Eabs,x=0.45,z=0.45,time=1 ../history.pe00000\*.nc@Eabs_ana,x=0.45,z=0.45,time=1 --nocont --aspect 1 --wsn 2 -sw:ifl=1
mv dcl_0001.png Eabs.png
gpview --exch --range 0:1.5 --overplot 2 ../history.pe00000\*.nc@Epot,x=0.45,z=0.45,time=1 ../history.pe00000\*.nc@Epot_ana,x=0.45,z=0.45,time=1 --nocont --aspect 1 --wsn 2 -sw:ifl=1
mv dcl_0001.png Epot.png
gpview --exch --overplot 2 ../history.pe00000\*.nc@CRGD_TOT,x=0.45,z=0.45,time=1 ../history.pe00000\*.nc@CRGD_TOT_ana,x=0.45,z=0.45,time=1 --nocont --aspect 1 --wsn 2 -sw:ifl=1
mv dcl_0001.png CRGD_TOT.png
