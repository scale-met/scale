#! /bin/bash -x
echo "==== POST PROCESS TO CREATE PICTURE ======="
.libs/profile_maker

mkdir figure
#---plot LWP by gnuplot
gnuplot << EOF
set terminal postscript color
set size square
set grid
set xlabel "Time [s]"
set ylabel "LWP [g/m2]"
plot "comp_time_evolv.txt" u 1:2 w l lw 3 t "LWP"
set output "LWP.ps"
replot
EOF
convert LWP.ps figure/LWP.jpeg
rm LWP.ps
#---plot ql by gnuplot
gnuplot << EOF
set terminal postscript color
set size square
set grid
set yrange [0:1200]
set xlabel "Liquid water mixing ratio [g/kg]"
set ylabel "Height [m]"
plot "comp_profile.txt" u 3:1 w l lw 5 t "LWP"
set output "ql.ps"
replot
EOF
convert ql.ps figure/ql.jpeg
rm ql.ps
#---plot qt by gnuplot
gnuplot << EOF
set terminal postscript color
set size square
set grid
set yrange [0:1200]
set xlabel "total water mixing ratio [g/kg]"
set ylabel "Height [m]"
plot "comp_profile.txt" u 2:1 w l lw 5 t "LWP"
set output "qt.ps"
replot
EOF
convert qt.ps figure/qt.jpeg
rm qt.ps

echo "draw QHYD at t=7200s"
gpview history.pe00000\*.nc@QHYD,time=7200 --mean x --nocont --aspect 1 --wsn 2
convert dcl.ps figure/qhyd_t7200.jpeg
rm dcl.ps
echo "==== POST PROCESS is finished ========="
