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

export OUTDIR=${HMDIR}/output/KH3/grd_1

mkdir -p ${OUTDIR}
cd ${OUTDIR}

########################################################################

# run
$BIN/$EXE GRID_DXYZ=40 GRID_KMAX=150 GRID_IMAX=25 GRID_JMAX=25 gridfile="../grid_40m_150x25x25" \
          PRC_nmax=72 PRC_NUM_X=18 PRC_NUM_Y=4 \
          step_str=1 step_end=25 \
          ../history
