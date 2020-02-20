set terminal png size 1280,640
set output "Mass_Ocean.png"

set xlabel "Step"
set ylabel "Mass[kg]"

set multiplot layout 2,1 title "Mass conservation check"

plot "mass_ocean.dat" using 1:2 ti "Mass (ICE)"  with linespoints pointsize 1.0, \
     "mass_ocean.dat" using 1:4 ti "Mass Change (ICE)" with linespoints pointsize 0.1
plot "mass_ocean.dat" using 1:3 ti "Mass Change (WATER)" with linespoints pointsize 0.2

unset multiplot