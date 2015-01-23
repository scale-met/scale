#!/bin/bash

cat << EOF > base.jobshell
#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################

##### RESOURCE SETTING
#-------------------------------------------------------------------------------
#PJM --rsc-list "rscgrp=${RSCGRP}"
#PJM --rsc-list "node=${TOTALNODE}"
#PJM --rsc-list "elapse=24:00:00"
#PJM --rsc-list "node-mem=13Gi"
#PJM --rsc-list "node-quota=29G"
#PJM --mpi assign-online-node
#PJM --mpi "use-rankdir"
#PJM --stg-transfiles all
#-------------------------------------------------------------------------------

##### COMMON SETTING
#-------------------------------------------------------------------------------
#PJM --stgin  "rank=* scale-les              %r:./"
#PJM --stgin  "rank=* run.d0*.conf           %r:./"
#PJM --stgin  "rank=* launch.conf            %r:./"
#PJM --stgin  "rank=* data/PARAG.29          %r:./"
#PJM --stgin  "rank=* data/PARAPC.29         %r:./"
#PJM --stgin  "rank=* data/VARDATA.RM29      %r:./"
#PJM --stgin  "rank=* data/cira.nc           %r:./"
#PJM --stgin  "rank=* data/MIPAS/*           %r:./MIPAS/"
#PJM --stgin  "rank=* data/param.bucket.conf %r:./"
#-------------------------------------------------------------------------------

##### FOR EACH DOMAINS
#-------------------------------------------------------------------------------
EOF

# make a process table
KEY=0
DOM=1
NPRC=${PRC_DOMAINS[`expr ${DOM} - 1`]}

i=0
while [ ${i} -lt ${TOTALNODE} ]
do
  FKEY=`printf "%06d" ${KEY}`

  echo "#PJM --stgin  \"rank="${i}" input/init_d0"${DOM}"_"${INITSEC}".pe"${FKEY}".nc %r:./\"" >> base.jobshell
  echo "#PJM --stgin  \"rank="${i}" input/topo_d0"${DOM}".pe"${FKEY}".nc %r:./\"" >> base.jobshell
  echo "#PJM --stgin  \"rank="${i}" input/landuse_d0"${DOM}".pe"${FKEY}".nc %r:./\"" >> base.jobshell
  if [ ${DOM} -eq 1 ]; then
    echo "#PJM --stgin  \"rank="${i}" input/boundary_d0"${DOM}".pe"${FKEY}".nc %r:./\"" >> base.jobshell
  fi

  KEY=`expr ${KEY} + 1`

  if [ ${KEY} -ge ${NPRC} ]; then
     KEY=0
     DOM=`expr ${DOM} + 1`
     NPRC=${PRC_DOMAINS[`expr ${DOM} - 1`]}
  fi

  i=`expr ${i} + 1`
done

cat << EOF >> base.jobshell
#-------------------------------------------------------------------------------

##### STAGE OUT
#-------------------------------------------------------------------------------
#PJM --stgout "rank=* %r:./* ./output/"
#PJM -j
#PJM -s
#-------------------------------------------------------------------------------

##### ENVIRONMENT SETTING
#-------------------------------------------------------------------------------
. /work/system/Env_base
export PARALLEL=8
export OMP_NUM_THREADS=8
#-------------------------------------------------------------------------------

##### RUN
mpiexec ./scale-les launch.conf || exit 1
EOF
