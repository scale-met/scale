#!/bin/sh
#-----------------------------------------
#   Test Script for XXX  (Ver.0.1)
#   2009/10/26 --- Ryuji Yoshida.
#-----------------------------------------

#rm -v *.txt

echo 'mpirun -n 1 ./driver run.d01.conf'
mpirun -n 1 ./driver run.d01.conf



#-----------------------------------------
#EOF
