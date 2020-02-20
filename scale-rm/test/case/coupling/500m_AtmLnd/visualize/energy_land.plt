set terminal png size 1280,640
set output "Energy_Land.png"

set xlabel "Step"
set ylabel "Energy[J]"

set multiplot layout 2,1 title "Energy conservation check"

plot "energy_land.dat" using 1:2 ti "Total internal Energy" with linespoints pointsize 1.0, \
     "energy_land.dat" using 1:6 ti "Energy change"         with linespoints pointsize 0.2

plot "energy_land.dat" using 1:6     ti "Energy change (Total)"  with linespoints pointsize 1.0, \
     "energy_land.dat" using 1:3     ti "Energy change (GH)"     with linespoints pointsize 0.2, \
     "energy_land.dat" using 1:4     ti "Energy change (ENGI)"   with linespoints pointsize 0.2, \
     "energy_land.dat" using 1:(-$5) ti "Energy change (Runoff)" with linespoints pointsize 0.2


unset multiplot
