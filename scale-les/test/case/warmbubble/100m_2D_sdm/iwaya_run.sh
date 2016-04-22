#! /bin/sh -x
#PBS -q SMALL
#PBS -o output_run.log
#PBS -e error_run.log
#PBS -l select=1:ncpus=4:mpiprocs=4
#PBS -N test_job
cd ${PBS_O_WORKDIR}
. /usr/share/modules/init/bash
module load mpt PrgEnv-intel hdf5  netcdf

#export FORT_FMT_RECL=400
# ulimit -s unlimited

# run
mpiexec_mpt dplace -s1 ./scale-les_init  init.conf  || exit
mpiexec_mpt dplace -s1 ./scale-les  run.conf  || exit

################################################################################
