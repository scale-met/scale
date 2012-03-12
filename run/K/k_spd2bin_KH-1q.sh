#! /bin/bash -x
#
# for K Computer
#
#PJM --rsc-list "node=1x1"
#PJM --rsc-list "elapse=01:00:00"
#PJM --rsc-list "node-mem=10Gi"
#PJM -s
#
. /work/system/Env_base_1.2.0-02
#
export PARALLEL=8
export OMP_NUM_THREADS=$PARALLEL
export LPG="lpgparm -s 32MB -d 32MB -h 32MB -t 32MB -p 32MB"
export fu08bf=1

export HMDIR=/work/scratch/user0171/scale3
export BIN=${HMDIR}/sbin/K
export EXE=spd2bin

export OUTDIR=${HMDIR}/output/MP_init

mkdir -p ${OUTDIR}
cd ${OUTDIR}

########################################################################

# run
$BIN/$EXE GRID_DXYZ=20 GRID_KMAX=336 GRID_IMAX=64 GRID_JMAX=64 gridfile="../init_warmbubble/grid_20m_336x64x64" \
          PRC_nmax=1 PRC_NUM_X=1 PRC_NUM_Y=1 \
          history
