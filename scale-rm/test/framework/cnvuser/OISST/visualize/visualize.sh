#! /bin/bash -x

### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

# boundary conditions
gpview sst.pe\*.nc@SST --mean time --nocont --aspect=1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf sst.png
rm -f dcl.pdf

# boundary conditions
gpview sst.pe\*.nc@SST --mean x,y --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf slice_sst.png
rm -f dcl.pdf
