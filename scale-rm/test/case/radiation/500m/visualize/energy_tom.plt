set terminal png size 1280,640
set output "Energy_TOM.png"

set xlabel "Step"
set ylabel "Energy[J]"

set multiplot layout 3,3 title "Energy conservation check (TOM)"

plot "energy_tom.dat" using 1:($3+$4) ti "Net TOM LW Flux"                    with linespoints pointsize 0.2
plot "energy_tom.dat" using 1:3       ti "TOM LW Upward   Flux(ENGTOM_LW_up)" with linespoints pointsize 0.2
plot "energy_tom.dat" using 1:4       ti "TOM LW Downward Flux(ENGTOM_LW_dn)" with linespoints pointsize 0.2
plot "energy_tom.dat" using 1:($5+$6) ti "Net TOM SW Flux"                    with linespoints pointsize 0.2
plot "energy_tom.dat" using 1:5       ti "TOM SW Upward   Flux(ENGTOM_SW_up)" with linespoints pointsize 0.2
plot "energy_tom.dat" using 1:6       ti "TOM SW Downward Flux(ENGTOM_SW_dn)" with linespoints pointsize 0.2
plot "energy_tom.dat" using 1:2       ti "Net TOM Radiation Flux(ENGTOM_RD)"  with linespoints pointsize 0.2
plot "energy_tom.dat" using 1:($3+$5) ti "TOM LW+SW Upward   Flux"            with linespoints pointsize 0.2
plot "energy_tom.dat" using 1:($4+$6) ti "TOM LW+SW Downward Flux"            with linespoints pointsize 0.2

unset multiplot
