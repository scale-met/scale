#! /bin/bash -x

### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

var__set=(ps u850)
rangeset=(auto auto)
time_set=(00000 01440 07200 11520 15840)

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[${i}]} == "auto" ]; then
      eddy=""
      range=""
   elif [ ${rangeset[${i}]} == "eddy" ]; then
      eddy="--eddy time"
      range="--range=-3:3 --eddy x"
   else
      eddy=""
      range="--range="${rangeset[${i}]}
   fi

   # snapshot
   for sec in ${time_set[@]}
   do
       gpview ${var}.nc@${var},time=${sec} --nocont ${range} --wsn 2 || exit
       convert -density 150 -rotate 90 +antialias dcl.pdf ${var}${sec}sec.png
       rm -f dcl.pdf
   done

   let i="${i} + 1"
done
