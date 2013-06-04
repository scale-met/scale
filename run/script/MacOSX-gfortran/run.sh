#!/bin/sh

bin=$1
conf=$2
pnum=$3

cmd="mpirun -n $pnum $bin $conf"


echo $cmd
$cmd
