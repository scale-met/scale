#! /bin/bash -x
#PJM -L "rscgrp=small"
#PJM -L "node=49"
#PJM -L "elapse=8:00:00"
#PJM -o output.file.txt
#PJM -e error.file.txt
#PJM --mail-list yousuke.sato@riken.jp
#PJM -m e 
#PJM -m b 
#PJM -j 

export PARALLEL=8
export OMP_NUM_THREADS=$PARALLEL
export fu08bf=10

mkdir /group2/gc24/c24045/dycom2_bin_8
cd /home/c24045/scale3/output/dycoms2_bin_8

mpiexec /home/c24045/scale3/bin/FX10/DYCOMS2RF01_uniso_grid_hbinw_aero_uniso276x8x8_hbinw_aero_ DYCOMS2RF01_uniso276x10x10_hbinw_.cnf

