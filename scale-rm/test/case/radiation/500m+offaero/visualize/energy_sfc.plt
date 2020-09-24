set terminal png size 1280,640
set output "Energy_SFC.png"

set xlabel "Step"
set ylabel "Energy[J]"

set multiplot layout 3,3 title "Energy conservation check (SFC)"

plot "energy_sfc.dat" using 1:($3+$4) ti "Net SFC LW Flux"                    with linespoints pointsize 0.2
plot "energy_sfc.dat" using 1:3       ti "SFC LW Upward   Flux(ENGSFC_LW_up)" with linespoints pointsize 0.2
plot "energy_sfc.dat" using 1:4       ti "SFC LW Downward Flux(ENGSFC_LW_dn)" with linespoints pointsize 0.2
plot "energy_sfc.dat" using 1:($5+$6) ti "Net SFC SW Flux"                    with linespoints pointsize 0.2
plot "energy_sfc.dat" using 1:5       ti "SFC SW Upward   Flux(ENGSFC_SW_up)" with linespoints pointsize 0.2
plot "energy_sfc.dat" using 1:6       ti "SFC SW Downward Flux(ENGSFC_SW_dn)" with linespoints pointsize 0.2
plot "energy_sfc.dat" using 1:2       ti "Net SFC Radiation Flux(ENGSFC_RD)"  with linespoints pointsize 0.2
plot "energy_sfc.dat" using 1:($3+$5) ti "SFC LW+SW Upward   Flux"            with linespoints pointsize 0.2
plot "energy_sfc.dat" using 1:($4+$6) ti "SFC LW+SW Downward Flux"            with linespoints pointsize 0.2

unset multiplot
