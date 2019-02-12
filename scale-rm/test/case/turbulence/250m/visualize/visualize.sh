#! /bin/bash -x

### Energy & Mass balance ###
echo "+visualize by gnuplot"
rm -f energy.dat energy_flx.dat mass.dat

while read -a line
do
   if [ ${line[0]} == "STEP=" ]; then
      echo ${line[1]} ${line[5]} ${line[6]}  ${line[7]}  ${line[8]}  ${line[9]}  >> energy.dat
      echo ${line[1]} ${line[9]} ${line[10]} ${line[11]} ${line[12]} ${line[13]} >> energy_flx.dat
      echo ${line[1]} ${line[3]} ${line[4]}                                      >> mass.dat
   fi
done < monitor.peall

gnuplot < ./visualize/energy.plt     || exit
gnuplot < ./visualize/energy_flx.plt || exit
gnuplot < ./visualize/mass.plt       || exit
rm -f energy.dat energy_flx.dat mass.dat



### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

var__set=(PT U W TKE_SMG)
rangeset=(eddy 4.7:5.3 -0.5:0.5 auto)
time_set=(00000 14400)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[${i}]} == "auto" ]; then
      range=""
   elif [ ${rangeset[${i}]} == "eddy" ]; then
      range="--eddy x"
   else
      range="--range="${rangeset[${i}]}
   fi

   # snapshot
   for sec in ${time_set[@]}
   do
       gpview history.pe\*.nc@${var},y=0,time=${sec} --nocont ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf ${var}${sec}sec.png
       rm -f dcl.pdf
       gpview history.pe\*.nc@${var},z=0,time=${sec} --nocont ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf sfc_${var}${sec}sec.png
       rm -f dcl.pdf
       gpview history.pe\*.nc@${var},z=1000,time=${sec} --nocont ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf z1000_${var}${sec}sec.png
       rm -f dcl.pdf
   done

   let i="${i} + 1"
done

var__set=(U V W W_PRIM2 W_PRIM3 PT PT_W_PRIM TKE_SMG TKE_RS Ri)
rangeset=(4.5:5.5 -0.05:0.05 -0.003:0.003 -0.1:1.2 -1:1 298:307 -30:100 -0.1:2.0 -0.1:1.5 -10:100)
time_set=(00000 07200 14400)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[${i}]} == "auto" ]; then
      range=""
   elif [ ${rangeset[${i}]} == "eddy" ]; then
      range="--eddy x"
   else
      range="--range="${rangeset[${i}]}
   fi

   # time series (profile)
   command=" --mean x,y --exch ${range} --wsn 2 --overplot 3 "
   for sec in ${time_set[@]}
   do
      command="${command} history.pe*.nc@${var},time=${sec} "
   done
   gpview $command
   convert -density 150 -rotate 90 +antialias dcl.pdf profile_${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done

