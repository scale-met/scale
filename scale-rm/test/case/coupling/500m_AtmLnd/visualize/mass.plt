set terminal png size 1280,640
set output "Mass.png"

set xlabel "Step"
set ylabel "Mass[kg]"

set multiplot layout 3,2 title "Mass conservation check"

plot "mass.dat" using 1:2              ti "Dry   Mass(QDRY)"                     with linespoints pointsize 0.2

plot "mass.dat" using 1:($3+$6+$7+$9)  ti "Water Mass(ATMOS+LAND)+RunOff"        with linespoints pointsize 1.0, \
     "mass.dat" using 1:($4-$5+$10+$9) ti "Water Mass change(ATMOS+LAND)+RunOff" with linespoints pointsize 0.2

plot "mass.dat" using 1:($3+$6+$7) ti "Water Mass(ATMOS+LAND)" with linespoints pointsize 1.0, \
     "mass.dat" using 1:3          ti "Water Mass(ATMOS)"      with linespoints pointsize 0.2, \
     "mass.dat" using 1:($6+$7)    ti "Water Mass(LAND)"       with linespoints pointsize 0.2

plot "mass.dat" using 1:($6+$7) ti "Land Mass(WATER+ICE)" with linespoints pointsize 1.0, \
     "mass.dat" using 1:6       ti "Land Mass(WATER)"     with linespoints pointsize 0.2,\
     "mass.dat" using 1:7       ti "Land Mass(ICE)"       with linespoints pointsize 0.2, \
     "mass.dat" using 1:10      ti "Land Mass change"      with linespoints pointsize 0.2

plot "mass.dat" using 1:3       ti "Water Mass(ATMOS)"        with linespoints pointsize 1.0, \
     "mass.dat" using 1:($4-$5) ti "Water Mass change(Total)" with linespoints pointsize 0.2, \
     "mass.dat" using 1:4       ti "Water Mass change(EVAP)"  with linespoints pointsize 0.2, \
     "mass.dat" using 1:(-$5)   ti "Water Mass change(PREC)"  with linespoints pointsize 0.2

plot "mass.dat" using 1:10      ti "Land Mass change(Total)"  with linespoints pointsize 1.0, \
     "mass.dat" using 1:8       ti "Land Mass change(Top)"    with linespoints pointsize 0.2, \
     "mass.dat" using 1:(-$9)   ti "Land Mass change(Runoff)" with linespoints pointsize 0.2


unset multiplot
