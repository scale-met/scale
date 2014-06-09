#! /bin/bash -x
#BSUB -q smp-small
#BSUB -n 1 
#BSUB -J sdm_test
#BSUB -o output_run.log
#BSUB -e error_run.log
#BSUB -P 1
#

#export FORT_FMT_RECL=400
# ulimit -s unlimited

# run
#mpirun -np ${LSB_DJOB_NUMPROC} dplace -s1 ./scale-les_init  init.conf  || exit
mpirun -np ${LSB_DJOB_NUMPROC} dplace -s1 ./scale-les  run.conf  || exit

################################################################################
