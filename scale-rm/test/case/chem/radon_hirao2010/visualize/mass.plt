set terminal png size 1280,640
set output "Mass.png"

set xlabel "Step"
set ylabel "Mass[kg] or [Bq]"

set multiplot layout 2,2 title "Mass conservation check"

plot "mass.dat" using 1:2 ti "Dry   Mass(QDRY)" with linespoints pointsize 0.2
plot "mass.dat" using 1:3 ti "Water Mass(QTOT)" with linespoints pointsize 0.2
plot "mass.dat" using 1:4 ti "Radon-222(RN222)" with linespoints pointsize 0.2

unset multiplot
