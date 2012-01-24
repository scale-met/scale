#!/bin/bash
i=1
while [ $i -le 7 ]
do
 fprofpx -p all -l0 -d pa$i -o output_prof_$i.csv -tcsv -Ihwm
 i=`expr $i + 1`
done
