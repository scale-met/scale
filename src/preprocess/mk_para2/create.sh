#!/bin/bash
FC="gfortran-mp-4.7" 
FFLAGS="-fconvert=big-endian -O3"

make
./mkpara

${FC} ${FFLAGS} -o make_inc_tracer make_tracer_inc.f90 
./make_inc_tracer

make clean
rm -f make_inc_tracer
