#! /bin/bash -x

### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

# boundary conditions
gpview landuse.pe\*.nc@FRAC_LAND  --nocont --aspect=4 --range=0:1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf frac_land.png
gpview landuse.pe\*.nc@FRAC_LAKE  --nocont --aspect=4 --range=0:1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf frac_lake.png
gpview landuse.pe\*.nc@FRAC_URBAN --nocont --aspect=4 --range=0:1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf frac_urban.png
gpview landuse.pe\*.nc@FRAC_PFT1  --nocont --aspect=4 --range=0:1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf frac_pft1.png
gpview landuse.pe\*.nc@FRAC_PFT2  --nocont --aspect=4 --range=0:1 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf frac_pft2.png
gpview landuse.pe\*.nc@INDEX_PFT1 --nocont --aspect=4 --levels 0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf index_pft1.png
gpview landuse.pe\*.nc@INDEX_PFT2 --nocont --aspect=4 --levels 0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf index_pft2.png
rm -f dcl.pdf
