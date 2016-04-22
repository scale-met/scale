#!/bin/sh
#-----------------------------------------
#   Test Script for XXX  (Ver.0.1)
#   2009/10/26 --- Ryuji Yoshida.
#-----------------------------------------

#rm -v *.txt

echo 'mpirun -n 13 ./launcher launch.conf'
mpirun -n 13 ./launcher launch.conf



#-----------------------------------------
#EOF
