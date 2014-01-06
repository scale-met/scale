#! /bin/bash -x
echo "==== POST PROCESS TO CREATE PICTURE ======="
make analy
.libs/profile_maker

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
convert LWP.ps -rotate 90 LWP.png
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
convert ql.ps -rotate 90 ql.png
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
convert qt.ps -rotate 90 qt.png
rm qt.ps
#---plot ww by gnuplot
gnuplot << EOF
set terminal postscript color
set size square
set grid
set yrange [0:1200]
set xlabel "Variance of w' [m2/s2]"
set ylabel "Height [m]"
plot "comp_profile.txt" u 6:1 w l lw 5 t "w'w'"
set output "ww.ps"
replot
EOF
convert ww.ps -rotate 90 ww.png
rm ww.ps
#---plot gr-TKE by gnuplot
gnuplot << EOF
set terminal postscript color
set size square
set grid
set yrange [0:1200]
set xlabel "Variance of w' [m2/s2]"
set ylabel "Height [m]"
plot "comp_profile.txt" u 4:1 w l lw 5 t "w'w'"
set output "gr_tke.ps"
replot
EOF
convert gr_tke.ps -rotate 90 gr_tke.png
rm gr_tke.ps


echo "draw QHYD at t=3600s"
gpview history.pe00000\*.nc@QHYD,time=3600 --mean x --nocont --aspect 1 --wsn 2
convert dcl.ps -rotate 90 qhyd_t3600.png
rm dcl.ps
echo "==== POST PROCESS is finished ========="
