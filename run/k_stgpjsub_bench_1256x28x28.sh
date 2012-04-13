#! /bin/bash -x
#
# for K computer
#
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=1x1"
#PJM --rsc-list "elapse=00:10:00"
#PJM --rsc-list "node-mem=10Gi"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin  "rank=* /data1/user0117/scale3/bin/K/scale3_init_1256x28x28_ndw6_     %r:./"
#PJM --stgin  "rank=* /data1/user0117/scale3/run/K/scale3_init_1256x28x28_ndw6_.cnf %r:./"
#PJM --stgin  "rank=* /data1/user0117/scale3/bin/K/scale3_1256x28x28_ndw6_          %r:./"
#PJM --stgin  "rank=* /data1/user0117/scale3/run/K/scale3_1256x28x28_ndw6_.cnf      %r:./"
#PJM --stgout "rank=* %r:./* /data1/user0117/scale3/output/scale3_1256x28x28_ndw6_/"
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=${PARALLEL}
export LPG="lpgparm -s 32MB -d 32MB -h 32MB -t 32MB -p 32MB"
export fprof="fpcoll -Ihwm,cpu -l0 -o Basic_Profile.txt -m 200000"
export fu08bf=1

# output directory
outdir=/data1/user0117/scale3/output/scale3_1256x28x28_ndw6_
mkdir -p ${outdir}

# run
         mpiexec ${LPG} ./scale3_init_1256x28x28_ndw6_ ./scale3_init_1256x28x28_ndw6_.cnf
${fprof} mpiexec ${LPG} ./scale3_1256x28x28_ndw6_      ./scale3_1256x28x28_ndw6_.cnf

# move to output directory
mv k_stgpjsub_bench_1256x28x28.sh.* ${outdir}
cp scale3_init_1256x28x28_ndw6_.cnf ${outdir}
cp scale3_1256x28x28_ndw6_.cnf      ${outdir}