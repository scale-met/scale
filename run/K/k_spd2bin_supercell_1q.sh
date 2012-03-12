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

export OUTDIR=${HMDIR}/output/supercell-1

mkdir -p ${OUTDIR}
cd ${OUTDIR}

########################################################################

# run
$BIN/$EXE GRID_DXYZ=400 GRID_KMAX=50 GRID_IMAX=50 GRID_JMAX=50 gridfile="../init_supercell/grid_400m_50x50x50" \
          PRC_nmax=64 PRC_NUM_X=8 PRC_NUM_Y=8 \
          history
