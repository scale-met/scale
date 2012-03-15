#! /bin/bash -x
#
# for 3 team (NST) cluster
#
#BSUB -q as
#BSUB -n 1
#BSUB -J spd2bin
#BSUB -o spd2bin_log
#BSUB -e spd2bin_error

export OMP_NUM_THREADS=1

export HMDIR=/home/yashiro/scale3
export BIN=${HMDIR}/sbin/NST-ifort
export EXE=spd2bin

export OUTDIR=${HMDIR}/output/supercell2/grd_1

mkdir -p ${OUTDIR}
cd ${OUTDIR}

########################################################################

# run
$BIN/$EXE GRID_DXYZ=400 GRID_KMAX=50 GRID_IMAX=50 GRID_JMAX=50 gridfile="../grid_400m_50x50x50" \
          PRC_nmax=64 PRC_NUM_X=8 PRC_NUM_Y=8 \
          step_str=1 step_end=30 \
          ../history
