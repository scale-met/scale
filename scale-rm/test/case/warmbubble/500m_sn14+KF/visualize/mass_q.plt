set terminal png size 1280,640
set output "Mass_q.png"

set xlabel "Step"
set ylabel "Mass[kg]"

set multiplot layout 3,3 title "Mass conservation check"

plot "mass_q.dat" using 1:3 ti "Water vapor(QV)"  with linespoints pointsize 0.2
plot "mass_q.dat" using 1:4 ti "Cloud Water(QC)"  with linespoints pointsize 0.2
plot "mass_q.dat" using 1:5 ti "Rain  Water(QR)"  with linespoints pointsize 0.2
plot "mass_q.dat" using 1:6 ti "Cloud Ice  (QI)"  with linespoints pointsize 0.2
plot "mass_q.dat" using 1:7 ti "Snow       (QS)"  with linespoints pointsize 0.2
plot "mass_q.dat" using 1:8 ti "Graupel    (QG)"  with linespoints pointsize 0.2
plot "mass_q.dat" using 1:2 ti "Water Mass(QTOT)" with linespoints pointsize 0.2

unset multiplot
