#! /bin/bash -x
#
# for K computer
#
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=1x1"
#PJM --rsc-list "elapse=00:30:00"
#PJM --rsc-list "node-mem=10Gi"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin  "rank=* /data1/user0117/scale3/bin/K/scale3_init_436x14x14_ndw6_     %r:./"
#PJM --stgin  "rank=* /data1/user0117/scale3/run/K/scale3_init_436x14x14_ndw6_.cnf %r:./"
#PJM --stgin  "rank=* /data1/user0117/scale3/bin/K/scale3_436x14x14_ndw6_          %r:./"
#PJM --stgin  "rank=* /data1/user0117/scale3/run/K/scale3_436x14x14_ndw6_.cnf      %r:./"
#PJM --stgout "rank=* %r:./*.* ./job%j/"
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=${PARALLEL}
export LPG="lpgparm -s 32MB -d 32MB -h 32MB -t 32MB -p 32MB"
export fprof="fpcoll -Ihwm,cpu -l0 -o Basic_Profile.txt -m 200000"
export fu08bf=1

# run
         mpiexec ${LPG} ./scale3_init_436x14x14_ndw6_ ./scale3_init_436x14x14_ndw6_.cnf
${fprof} mpiexec ${LPG} ./scale3_436x14x14_ndw6_      ./scale3_436x14x14_ndw6_.cnf
