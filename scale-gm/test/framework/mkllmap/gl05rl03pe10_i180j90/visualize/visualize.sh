#! /bin/bash -x

### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

gpview sample1.nc@sample1 --nocont --itr=31 --map=coast_world --map_axis=0,90,0  --tone=b --int=0.5 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf rank_np.png
rm -f dcl.pdf

gpview sample1.nc@sample1 --nocont --itr=31 --map=coast_world --map_axis=0,-90,0 --tone=b --int=0.5 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf rank_sp.png
rm -f dcl.pdf

gpview sample2.nc@sample2 --nocont --itr=31 --map=coast_world --map_axis=0,90,0  --tone=b           --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf region_np.png
rm -f dcl.pdf

gpview sample2.nc@sample2 --nocont --itr=31 --map=coast_world --map_axis=0,-90,0 --tone=b           --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf region_sp.png
rm -f dcl.pdf

gpview sample3.nc@sample3 --nocont --itr=31 --map=coast_world --map_axis=0,90,0  --tone=b --int=0.5 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf gridi_np.png
rm -f dcl.pdf

gpview sample3.nc@sample3 --nocont --itr=31 --map=coast_world --map_axis=0,-90,0 --tone=b --int=0.5 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf gridi_sp.png
rm -f dcl.pdf

gpview sample4.nc@sample4 --nocont --itr=31 --map=coast_world --map_axis=0,90,0  --tone=b --int=0.5 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf gridj_np.png
rm -f dcl.pdf

gpview sample4.nc@sample4 --nocont --itr=31 --map=coast_world --map_axis=0,-90,0 --tone=b --int=0.5 --wsn 2 || exit
convert -density 150 -rotate 90 +antialias dcl.pdf gridj_sp.png
rm -f dcl.pdf
