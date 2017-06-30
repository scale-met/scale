set terminal png size 1280,640
set output "Energy_TOA.png"

set xlabel "Step"
set ylabel "Energy[J]"

set multiplot layout 3,3 title "Energy conservation check (TOA)"

plot "energy_toa.dat" using 1:($3+$4) ti "Net TOA LW Flux"                    with linespoints pointsize 0.2
plot "energy_toa.dat" using 1:3       ti "TOA LW Upward   Flux(ENGTOA_LW_up)" with linespoints pointsize 0.2
plot "energy_toa.dat" using 1:4       ti "TOA LW Downward Flux(ENGTOA_LW_dn)" with linespoints pointsize 0.2
plot "energy_toa.dat" using 1:($5+$6) ti "Net TOA SW Flux"                    with linespoints pointsize 0.2
plot "energy_toa.dat" using 1:5       ti "TOA SW Upward   Flux(ENGTOA_SW_up)" with linespoints pointsize 0.2
plot "energy_toa.dat" using 1:6       ti "TOA SW Downward Flux(ENGTOA_SW_dn)" with linespoints pointsize 0.2
plot "energy_toa.dat" using 1:2       ti "Net TOA Radiation Flux(ENGTOA_RD)"  with linespoints pointsize 0.2
plot "energy_toa.dat" using 1:($3+$5) ti "TOA LW+SW Upward   Flux"            with linespoints pointsize 0.2
plot "energy_toa.dat" using 1:($4+$6) ti "TOA LW+SW Downward Flux"            with linespoints pointsize 0.2

unset multiplot
