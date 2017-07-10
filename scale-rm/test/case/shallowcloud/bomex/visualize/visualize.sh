#! /bin/bash -x

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
done < monitor.pe000000

gnuplot < ./visualize/energy.plt     || exit
gnuplot < ./visualize/energy_flx.plt || exit
gnuplot < ./visualize/mass.plt       || exit
gnuplot < ./visualize/mass_q.plt     || exit
rm -f energy.dat energy_flx.dat mass.dat mass_q.dat



### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

var__set=(LWP)
rangeset=(auto)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time series
   gpview history.pe\*.nc@${var},time=7500:21600 --mean x,y --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf ${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done

var__set=(TKE_RS TKE_SMG)
rangeset=(auto)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time series
   gpview history.pe\*.nc@${var},time=7500:21600 --mean z,x,y --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf seq_${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done

var__set=(QTOT TKE_RS TKE_SMG PT   U    V    W_PRIM2 W_PRIM3 PT_W_PRIM)
rangeset=(auto auto   auto    auto auto auto auto    auto    auto     )

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time average
       gpview history.pe\*.nc@${var},time=18000:21600 --mean x,y,time --exch ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf prof_${var}.png
       rm -f dcl.pdf

   let i="${i} + 1"
done
