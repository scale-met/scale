#! /bin/bash -x

### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

# boundary conditions
gpview metrics.pe\*.nc@lat       --nocont --aspect=1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf lat.png
gpview metrics.pe\*.nc@lon       --nocont --aspect=1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf lon.png
gpview metrics.pe\*.nc@MAPF_X_XY --nocont --aspect=1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf MAPF_X_XY.png
gpview metrics.pe\*.nc@MAPF_Y_XY --nocont --aspect=1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf MAPF_Y_XY.png
gpview metrics.pe\*.nc@ROTC_COS  --nocont --aspect=1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf ROTC_COS.png
gpview metrics.pe\*.nc@ROTC_SIN  --nocont --aspect=1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf ROTC_SIN.png
gpview metrics.pe\*.nc@X_XY      --nocont --aspect=1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf x_xy.png
gpview metrics.pe\*.nc@Y_XY      --nocont --aspect=1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf y_xy.png
gpview metrics.pe\*.nc@distance  --nocont --aspect=1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf distance.png
rm -f dcl.pdf
