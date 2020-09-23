#! /bin/bash -x

### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

for domain in d01 d02
do

# metrics
gpview topo_${domain}.pe\*.nc@lat  --nocont --aspect=1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf lat_${domain}.png
gpview topo_${domain}.pe\*.nc@lon  --nocont --aspect=1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf lon_${domain}.png

# boundary conditions
gpview topo_${domain}.pe\*.nc@topo --nocont --aspect=1 --wsn 2 --srange=10: || exit
convert -density 150 -rotate 90 +antialias dcl.pdf topo_${domain}.png

gpview landuse_${domain}.pe\*.nc@FRAC_LAND  --nocont --aspect=1 --range=0:1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf frac_land_${domain}.png
# gpview landuse_${domain}.pe\*.nc@FRAC_LAKE  --nocont --aspect=1 --range=0:1 --wsn 2 || exit
# convert -density 150 -rotate 90 +antialias dcl.pdf frac_lake_${domain}.png
gpview landuse_${domain}.pe\*.nc@FRAC_URBAN --nocont --aspect=1 --range=0:1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf frac_urban_${domain}.png
# gpview landuse_${domain}.pe\*.nc@FRAC_PFT1  --nocont --aspect=1 --range=0:1 --wsn 2 || exit
# convert -density 150 -rotate 90 +antialias dcl.pdf frac_pft1_${domain}.png
# gpview landuse_${domain}.pe\*.nc@FRAC_PFT2  --nocont --aspect=1 --range=0:1 --wsn 2 || exit
# convert -density 150 -rotate 90 +antialias dcl.pdf frac_pft2_${domain}.png
# gpview landuse_${domain}.pe\*.nc@INDEX_PFT1 --nocont --aspect=1 --levels 0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5 --wsn 2 || exit
# convert -density 150 -rotate 90 +antialias dcl.pdf index_pft1_${domain}.png
# gpview landuse_${domain}.pe\*.nc@INDEX_PFT2 --nocont --aspect=1 --levels 0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5 --wsn 2 || exit
# convert -density 150 -rotate 90 +antialias dcl.pdf index_pft2_${domain}.png
rm -f dcl.pdf

done
