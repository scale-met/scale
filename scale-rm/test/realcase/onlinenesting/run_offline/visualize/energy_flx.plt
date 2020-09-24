set terminal png size 1280,640
set output "Energy_FLX.png"

set xlabel "Step"
set ylabel "Energy[J]"

set multiplot layout 3,2 title "Flux conservation check"

plot "energy_flx.dat" using 1:2          ti "Total Flux(ENGFLXT)"                 with linespoints pointsize 0.2
plot "energy_flx.dat" using 1:($3+$4+$5) ti "Net SFC Flux"                        with linespoints pointsize 0.2
plot "energy_flx.dat" using 1:3          ti "Sensible heat Flux (ENGFLX\\_SH)"    with linespoints pointsize 0.2
plot "energy_flx.dat" using 1:4          ti "Latent   heat Flux (ENGFLX\\_LH)"    with linespoints pointsize 0.2
plot "energy_flx.dat" using 1:5          ti "Net SFC Radiation Flux(ENGSFC\\_RD)" with linespoints pointsize 0.2
plot "energy_flx.dat" using 1:6          ti "Net TOM Radiation Flux(ENGTOM\\_RD)" with linespoints pointsize 0.2

unset multiplot
