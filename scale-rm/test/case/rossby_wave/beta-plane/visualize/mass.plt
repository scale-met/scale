set terminal png size 1280,640
set output "Mass.png"

set xlabel "Step"
set ylabel "Mass[kg]"

set multiplot layout 2,2 title "Mass conservation check"

plot "mass.dat" using 1:2 ti "Dry   Mass(QDRY)"   with linespoints pointsize 0.2

unset multiplot