set terminal png size 1280,640
set output "Energy.png"

set xlabel "Step"
set ylabel "Energy[J]"

set multiplot layout 2,2 title "Energy conservation check"

plot "energy.dat" using 1:3  ti "Net SFC Radiation Flux(ENGSFC\\_RD)" with linespoints pointsize 0.2
plot "energy.dat" using 1:4  ti "Net TOM Radiation Flux(ENGTOM\\_RD)" with linespoints pointsize 0.2
plot "energy.dat" using 1:2  ti "Total Radiation Flux(ENGFLXT)"       with linespoints pointsize 0.2

unset multiplot
