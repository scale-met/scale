set terminal png size 1280,640
set output "Energy_FLX.png"

set xlabel "Step"
set ylabel "Energy[J]"

set multiplot layout 2,2 title "Flux conservation check"

plot "energy_flx.dat" using 1:2 ti "Total Flux(ENGFLXT)" with linespoints pointsize 0.2

plot "energy_flx.dat" using 1:($3+$4) ti "Net SFC heat Flux"      with linespoints pointsize 0.2, \
     "energy_flx.dat" using 1:($5-$6) ti "Net SFC ENGI Flux"      with linespoints pointsize 0.2, \
     "energy_flx.dat" using 1:7       ti "Net SFC Radiation Flux" with linespoints pointsize 0.2, \
     "energy_flx.dat" using 1:8       ti "Net TOA Radiation Flux" with linespoints pointsize 0.2

plot "energy_flx.dat" using 1:3                ti "Sensible heat Flux (ENGFLX\\_SH)"    with linespoints pointsize 0.2, \
     "energy_flx.dat" using 1:4                ti "Latent   heat Flux (ENGFLX\\_LH)"    with linespoints pointsize 0.2

plot "energy_flx.dat" using 1:5                ti "Evapolation ENGI Flux (ENGFLX\\_EVAP)"  with linespoints pointsize 0.2, \
     "energy_flx.dat" using 1:(-$6)            ti "Precipitation ENGI Flux (ENGFLX\\_PREC)"  with linespoints pointsize 0.2

unset multiplot
