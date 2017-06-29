#!/bin/bash -x

cp init.conf init.orig.conf
rm ./*.png

for LAT in 0N 90N 90S
do
for LON in 0E 90E 90W 450E 450W
do
   if [ -f init_${LAT}-${LON}.conf ]; then
      cp -v init_${LAT}-${LON}.conf init.conf

      make run || continue

      bash eval.sh || continue

      mv lat.png       lat_${LAT}-${LON}.png
      mv lon.png       lon_${LAT}-${LON}.png
      mv MAPF_X_XY.png MAPF_X_XY_${LAT}-${LON}.png
      mv MAPF_Y_XY.png MAPF_Y_XY_${LAT}-${LON}.png
      mv ROTC_COS.png  ROTC_COS_${LAT}-${LON}.png
      mv ROTC_SIN.png  ROTC_SIN_${LAT}-${LON}.png
      mv x_xy.png      x_xy_${LAT}-${LON}.png
      mv y_xy.png      y_xy_${LAT}-${LON}.png
      mv distance.png  distance_${LAT}-${LON}.png
   fi
done
done
