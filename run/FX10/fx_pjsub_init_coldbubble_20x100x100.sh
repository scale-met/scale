#! /bin/bash -x
#
# for OAKLEAF-FX
#
#PJM -L "rscgrp=short"
#PJM -L "node=3x2"
#PJM -L "elapse=00:10:00"
#PJM -j

export PARALLEL=8
export OMP_NUM_THREADS=$PARALLEL
export fu08bf=1

export HMDIR=/home/c18001/scale3
export BIN=${HMDIR}/bin/FX10
export EXE=init_coldbubble_20x100x100_ndw6_

export OUTDIR=${HMDIR}/output/init_coldbubble_1km

mkdir -p ${OUTDIR}
cd ${OUTDIR}
rm -rf ./prof

########################################################################
cat << End_of_SYSIN > ${OUTDIR}/${EXE}.cnf

#####
#
# Scale3 init_coldbubble configulation
#
#####

&PARAM_PRC
 PRC_NUM_X       = 3,
 PRC_NUM_Y       = 2,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_GRID
 GRID_OUT_BASENAME = "grid_1000m_20x100x100",
/

&PARAM_GEOMETRICS
 GEOMETRICS_startlonlat  = 120.D0, 30.D0,
 GEOMETRICS_rotation     = 10.D0,
 GEOMETRICS_OUT_BASENAME = "",
/

&PARAM_COMM
 COMM_total_doreport  = .true.,
 COMM_total_globalsum = .true.,
/

&PARAM_ATMOS_VARS
 ATMOS_RESTART_OUTPUT         = .true.,
 ATMOS_RESTART_OUT_BASENAME   = "init_coldbubble",
/

&PARAM_MKEXP_COLDBUBBLE
 ZC_BBL =   5.D3,
 XC_BBL = 100.D3,
 YC_BBL = 100.D3,
 ZR_BBL =   2.D3,
 XR_BBL =   2.D3,
 YR_BBL =   2.D3,
/

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
fipp -C -Ihwm -d prof mpiexec $BIN/$EXE $EXE.cnf
echo "job ${RUNNAME} end     at " `date`

exit
