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
#PJM --stgin  "rank=* /data1/user0117/scale3/run/K/1x1/scale3_init1_005h_hydrostatic.cnf %r:./"
#PJM --stgin  "rank=* /data1/user0117/scale3/run/K/1x1/scale3_test1_005h_hydrostatic.cnf %r:./"
#PJM --stgin  "rank=* /data1/user0117/scale3/bin/K/scale3_init_1256x32x32_ndw6_          %r:./"
#PJM --stgin  "rank=* /data1/user0117/scale3/bin/K/scale3_1256x32x32_ndw6_               %r:./"
#PJM --stgout "rank=* %r:./* /data1/user0117/scale3/output/scale3_test1_005h_hydrostatic_1x1/"
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=${PARALLEL}
export LPG="lpgparm -s 32MB -d 32MB -h 32MB -t 32MB -p 32MB"
export fprof="fpcoll -Ihwm,cpu -l0 -o Basic_Profile.txt -m 200000"
export fu08bf=1

outdir=/data1/user0117/scale3/output/scale3_test1_005h_hydrostatic_1x1
mkdir -p ${outdir}

# run
         mpiexec ${LPG} ./scale3_init_1256x32x32_ndw6_ ./scale3_init1_005h_hydrostatic.cnf
${fprof} mpiexec ${LPG} ./scale3_1256x32x32_ndw6_      ./scale3_test1_005h_hydrostatic.cnf

mv k_stgpjsub_test1_005h_hydrostatic.sh.* ${outdir}
cp scale3_init1_005h_hydrostatic.cnf      ${outdir}
cp scale3_test1_005h_hydrostatic.cnf      ${outdir}
