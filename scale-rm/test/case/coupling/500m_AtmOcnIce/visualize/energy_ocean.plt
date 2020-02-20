set terminal png size 1280,640
set output "Energy_Ocean.png"

set xlabel "Step"
set ylabel "Energy[J]"

set multiplot layout 2,2 title "Energy conservation check"

plot "energy_ocean.dat" using 1:2 ti "ENG (Ocean)"           with linespoints pointsize 1.0, \
     "energy_ocean.dat" using 1:9 ti "ENG change (Ocean)"    with linespoints pointsize 0.1
plot "energy_ocean.dat" using 1:3 ti "ENG (Ice)"             with linespoints pointsize 1.0, \
     "energy_ocean.dat" using 1:10 ti "ENG change (Ice)"     with linespoints pointsize 0.1
plot "energy_ocean.dat" using 1:9 ti "ENG change (Ocean)"    with linespoints pointsize 0.2, \
     "energy_ocean.dat" using 1:6 ti "ENG change (GH)"       with linespoints pointsize 0.2, \
     "energy_ocean.dat" using 1:7 ti "ENG change (ENGI)"     with linespoints pointsize 0.2, \
     "energy_ocean.dat" using 1:8 ti "ENG change (SUPL)"     with linespoints pointsize 0.2
plot "energy_ocean.dat" using 1:10 ti "ENG change (Ice)"     with linespoints pointsize 0.2, \
     "energy_ocean.dat" using 1:4 ti "ENG change (GH top)"   with linespoints pointsize 0.2, \
     "energy_ocean.dat" using 1:5 ti "ENG change (ENGI top)" with linespoints pointsize 0.2, \
     "energy_ocean.dat" using 1:(-$6) ti "ENG change (GH bottom)"   with linespoints pointsize 0.2, \
     "energy_ocean.dat" using 1:(-$7) ti "ENG change (ENGI bottom)" with linespoints pointsize 0.2, \

unset multiplot
