#! /bin/bash -x
#
# for 3 team (NST) cluster
#
#BSUB -q bs
#BSUB -n 1
#BSUB -J spd2bin
#BSUB -o spd2bin_log
#BSUB -e spd2bin_error

export OMP_NUM_THREADS=1

export HMDIR=/home/yashiro/scale3
export BIN=${HMDIR}/sbin/NST-ifort
export EXE=spd2bin

export OUTDIR=${HMDIR}/output/supercell

cd ${OUTDIR}

########################################################################

# run
$BIN/$EXE GRID_DXYZ=400 GRID_KMAX=50 GRID_IMAX=100 GRID_JMAX=100 gridfile="grid_400m_50x100x100" \
          PRC_nmax=16 PRC_NUM_X=4 PRC_NUM_Y=4 \
          step_end=20 \
          history
