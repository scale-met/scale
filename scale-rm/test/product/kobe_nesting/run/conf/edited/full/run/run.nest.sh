#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################

##### RESOURCE SETTING
#-------------------------------------------------------------------------------
#PJM --rsc-list "rscgrp=large"
#PJM --rsc-list "node=15112"
#PJM --rsc-list "elapse=01:00:00"
#PJM --mpi      "shape=36"
#PJM --mpi      "proc=36"
#PJM --mpi      "use-rankdir"
#PJM --stg-transfiles all
#-------------------------------------------------------------------------------

##### FOR TOP PARENT DOMAIN
#-------------------------------------------------------------------------------
#PJM --stgin  "rank=* /data/ra000006/a00000/scale/bin/scale-rm        %r:./"
#PJM --stgin  "rank=* ../config/run.d01.conf                               %r:./"
#PJM --stgin  "rank=* ../init/domain_01/init_d01_00023587200.000.pe%06r.nc %r:./"
#PJM --stgin  "rank=* ../init/domain_01/boundary_d01.pe%06r.nc             %r:./"
#PJM --stgin  "rank=* ../input/domain_01/topo_d01.pe%06r.nc                %r:./"
#PJM --stgin  "rank=* ../input/domain_01/landuse_d01.pe%06r.nc             %r:./"
#PJM --stgin  "rank=* /data/ra000006/SCALE/database/rad/PARAG.29           %r:./"
#PJM --stgin  "rank=* /data/ra000006/SCALE/database/rad/PARAPC.29          %r:./"
#PJM --stgin  "rank=* /data/ra000006/SCALE/database/rad/VARDATA.RM29       %r:./"
#PJM --stgin  "rank=* /data/ra000006/SCALE/database/rad/cira.nc            %r:./"
#PJM --stgin  "rank=* /data/ra000006/SCALE/database/rad/MIPAS/*      %r:./MIPAS/"
#-------------------------------------------------------------------------------

##### FOR PARENT DOMAIN
#-------------------------------------------------------------------------------
#>>> DOMAIN 2
#PJM --stgin  "rank=2 ../config/run.d02.conf                             %r:../"
#PJM --stgin  "rank=2 ../init/domain_02/init_d02_00023587200.000.pe*.nc  %r:../"
#--- --stgin  "rank=2 ../init/domain_02/boundary_d02.pe*.nc              %r:../"
#PJM --stgin  "rank=2 ../input/domain_02/topo_d02.pe*.nc                 %r:../"
#PJM --stgin  "rank=2 ../input/domain_02/landuse_d02.pe*.nc              %r:../"
#>>> DOMAIN 3
#PJM --stgin  "rank=3 ../config/run.d03.conf                             %r:../"
#PJM --stgin  "rank=3 ../init/domain_03/init_d03_00023587200.000.pe*.nc  %r:../"
#--- --stgin  "rank=3 ../init/domain_03/boundary_d03.pe*.nc              %r:../"
#PJM --stgin  "rank=3 ../input/domain_03/topo_d03.pe*.nc                 %r:../"
#PJM --stgin  "rank=3 ../input/domain_03/landuse_d03.pe*.nc              %r:../"
#>>> DOMAIN 4
#PJM --stgin  "rank=4 ../config/run.d04.conf                             %r:../"
#PJM --stgin  "rank=4 ../init/domain_04/init_d04_00023587200.000.pe*.nc  %r:../"
#--- --stgin  "rank=4 ../init/domain_04/boundary_d04.pe*.nc              %r:../"
#PJM --stgin  "rank=4 ../input/domain_04/topo_d04.pe*.nc                 %r:../"
#PJM --stgin  "rank=4 ../input/domain_04/landuse_d04.pe*.nc              %r:../"

#>>> COMMON FILES
#PJM --stgin  "rank=0 /data/ra000006/a00000/scale/bin/scale-rm       0:../"
#PJM --stgin  "rank=0 /data/ra000006/SCALE/database/rad/PARAG.29          0:../"
#PJM --stgin  "rank=0 /data/ra000006/SCALE/database/rad/PARAPC.29         0:../"
#PJM --stgin  "rank=0 /data/ra000006/SCALE/database/rad/VARDATA.RM29      0:../"
#PJM --stgin  "rank=0 /data/ra000006/SCALE/database/rad/cira.nc           0:../"
#PJM --stgin  "rank=0 /data/ra000006/SCALE/database/rad/MIPAS/*     0:../MIPAS/"
#-------------------------------------------------------------------------------

##### STAGE OUT
#-------------------------------------------------------------------------------
#PJM --stgout "rank=* %r:./*       ./output/"
#PJM --stgout "rank=* %r:../*      ./output/"
#PJM --stgout "rank=* %r:./prof/*  ./output/prof/"
#PJM -j
#PJM -s
#-------------------------------------------------------------------------------

##### ENVIRONMENT SETTING
#-------------------------------------------------------------------------------
. /work/system/Env_base
export PARALLEL=8
export OMP_NUM_THREADS=8

fprof="fipp -C -Srange -Icall,hwm -d prof"
#-------------------------------------------------------------------------------


##### RUN
#ls -lha ./
#ls -lha ../

${fprof} mpiexec -n 36 ./scale-rm run.d01.conf || exit 1

