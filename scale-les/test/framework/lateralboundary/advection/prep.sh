#!/bin/sh
#-----------------------------------------
#   Preprocess Running Script for Boundary Experiments (Ver.0.1)
#   2014/06/09 --- Ryuji Yoshida.
#-----------------------------------------
home=`pwd`

cd ../prep_boundary_exp/
make
cd $home

ln -s ../prep_boundary_exp/bound_exp ./
./bound_exp

#-----------------------------------------
#EOF
