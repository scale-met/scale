#! /bin/bash -x

### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

var__set=(PT_urb QA_urb UA_urb SWD_urb RAIN_urb URBAN_SFC_TEMP URBAN_TC URBAN_T2 URBAN_QC URBAN_Q2 URBAN_UC URBAN_U10 URBAN_V10 URBAN_RAINR URBAN_RNgrd URBAN_ROFF)
rangeset=(290:320 0.01:0.02 0:5 0:1000 0:0.1 290:320 290:320 290:320 0:0.1 0:0.1 0:5 0:5 0:5 0:0.3 -200:800 0:5)
time_set=

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time series
   gpview history.pe000000.nc@${var}, --mean x,y ${range} --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf ${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done

var__set=(URBAN_TRL URBAN_TBL URBAN_TGL)
rangeset=(290:320 290:320 290:320)
time_set=

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time series
   gpview history.pe000000.nc@${var} --mean x,y ${range} --nocont --exch --wsn 2 || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf ${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done

var__set=(FLX_total TEMP   SHFLX   LHFLX   GHFLX   NetRad )
titl_set=(LH,SH,GH TR,TB,TG SHR,SHB,SHG LHR,LHB,LHG GHR,GHB,GHG RNR,RNB,RNG)
var1_set=(URBAN_SFLX_LH URBAN_TR URBAN_SHR URBAN_LHR URBAN_GHR URBAN_RNR)
var2_set=(URBAN_SFLX_SH URBAN_TB URBAN_SHB URBAN_LHB URBAN_GHB URBAN_RNB)
var3_set=(URBAN_SFLX_GH URBAN_TG URBAN_SHG URBAN_LHG URBAN_GHG URBAN_RNG)
rangeset=(-400:400 290:320 -400:400 -400:400 -400:400 -200:800)
time_set=

i=0
for var in ${var__set[@]}
do
   if [ ${rangeset[$i]} == "auto" ]; then
      range=""
   else
      range="--range="${rangeset[$i]}
   fi

   # time series
   gpview --wsn 2 --title ${titl_set[$i]} --mean x,y ${range} --overplot 3 history.pe000000.nc@${var1_set[$i]} \
                                                                           history.pe000000.nc@${var2_set[$i]} \
                                                                           history.pe000000.nc@${var3_set[$i]} || exit
   convert -density 150 -rotate 90 +antialias dcl.pdf ${var}.png
   rm -f dcl.pdf

   let i="${i} + 1"
done
