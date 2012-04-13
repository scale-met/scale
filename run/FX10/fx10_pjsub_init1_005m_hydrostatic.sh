#! /bin/bash -x
#
# for OAKLEAF-FX
#
#PJM -L "rscgrp=short"
#PJM -L "node=1"
#PJM -L "elapse=00:10:00"
#PJM -j
#PJM -s
# export postfix=${1:-20120413}

export PARALLEL=8
export OMP_NUM_THREADS=$PARALLEL
export fu08bf=10

export HMDIR=/home/c18001/scale3
export BIN=${HMDIR}/bin/FX10
export EXE=scale3_init_436x14x14_ndw6_${postfix}

export OUTDIR=${HMDIR}/output/${postfix}/init1_005m_hydrostatic

mkdir -p ${OUTDIR}
cd ${OUTDIR}

########################################################################
cat << End_of_SYSIN > ${OUTDIR}/${EXE}.cnf

#####
#
# Scale3 mkinit configulation
#
#####

&PARAM_PRC
 PRC_NUM_X       = 1,
 PRC_NUM_Y       = 1,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_GRID
 GRID_OUT_BASENAME = "grid_005m_436x14x14",
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
 ATMOS_RESTART_OUT_BASENAME   = "init1_hydrostatic",
/

&PARAM_MKINIT
 MKINIT_initname = "PLANESTATE",
/

&PARAM_MKINIT_PLANESTATE
 SFC_RH       = 50.D0,
 ENV_RH       = 50.D0,
 ENV_U        = 10.D0,
 ENV_V        = 10.D0,
 RANDOM_THETA =  1.D-1,
/

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
mpiexec $BIN/$EXE $EXE.cnf
echo "job ${RUNNAME} end     at " `date`

exit
