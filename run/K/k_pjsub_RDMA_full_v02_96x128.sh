#! /bin/bash -x
#
# for K Computer
#
#PJM --rsc-list "rscgrp=large"
#PJM --rsc-list "node=128x96"
#PJM --rsc-list "elapse=00:50:00"
#PJM --rsc-list "node-mem=10Gi"
#PJM -s
#
. /work/system/Env_base_1.2.0-02
#
export PARALLEL=8
export OMP_NUM_THREADS=$PARALLEL
export LPG="lpgparm -s 32MB -d 32MB -h 32MB -t 32MB -p 32MB"
export fu08bf=10

export HMDIR=/work/scratch/user0171/scale3
export BIN=${HMDIR}/bin/K
export EXE=RDMA_full_v02

export OUTDIR=${HMDIR}/output/RDMA_full_v02_128x96

mkdir -p ${OUTDIR}
cd ${OUTDIR}

########################################################################
cat << End_of_SYSIN > ${OUTDIR}/${EXE}.cnf

#####
#
# Scale3 benchmark configuration
#
#####

&PARAM_PRC
 PRC_NUM_X       = 128,
 PRC_NUM_Y       = 96,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_COMM
 COMM_vsize_max = 40,
/

&PARAM_TIME
 TIME_STARTDATE             = 2000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
 TIME_DURATION              = 3.0D0,
 TIME_DURATION_UNIT         = "SEC",
 TIME_DT                    = 0.6D0,
 TIME_DT_UNIT               = "SEC",
 TIME_DT_ATMOS_DYN          = 0.03D0,
 TIME_DT_ATMOS_DYN_UNIT     = "SEC",
 TIME_NSTEP_ATMOS_DYN       = 20,
 TIME_DT_ATMOS_PHY_MP       = 0.6D0,
 TIME_DT_ATMOS_PHY_MP_UNIT  = "SEC",
 TIME_DT_ATMOS_RESTART      = 3.0D0,
 TIME_DT_ATMOS_RESTART_UNIT = "SEC",
/

&PARAM_GRID
 GRID_OUT_BASENAME = '',
 GRID_DXYZ         = 20.D0,
 GRID_KMAX         = 336,
 GRID_IMAX         = 64,
 GRID_JMAX         = 64,
 GRID_BUFFER_DZ    = 6.0D3,
 GRID_BUFFFACT     = 1.1D0,
/

&PARAM_ATMOS
 ATMOS_TYPE_DYN    = "fent_fct",
 ATMOS_TYPE_PHY_TB = "smagorinsky",
 ATMOS_TYPE_PHY_MP = "NDW6",
 ATMOS_TYPE_PHY_RD = "mstrnX",
/

&PARAM_ATMOS_VARS
 ATMOS_QTRC_NMAX              = 11,
 ATMOS_RESTART_IN_BASENAME    = "${HMDIR}/output/init_warmbubble/init_warmbubble_63072000000.000",
 ATMOS_RESTART_OUTPUT         = .false.,
 ATMOS_RESTART_OUT_BASENAME   = "check_warmbubble",
 ATMOS_RESTART_CHECK          = .false.,
 ATMOS_RESTART_CHECK_BASENAME = "${HMDIR}/data/warmbubble/check_warmbubble_63072000003.000",
/

&PARAM_ATMOS_REFSTATE
 ATMOS_REFSTATE_TYPE         = "UNIFORM",
 ATMOS_REFSTATE_POTT_UNIFORM = 300.D0
/

&PARAM_ATMOS_BOUNDARY
 ATMOS_BOUNDARY_OUT_BASENAME = "boundary",
 ATMOS_BOUNDARY_VALUE_VELX   =  0.D0,
 ATMOS_BOUNDARY_TAUZ         = 10.D0,
/

&PARAM_ATMOS_DYN
 ATMOS_DYN_NUMERICAL_DIFF = 1.D-3,
/

&PARAM_OCEAN
 OCEAN_TYPE = "FIXEDSST",
/

&PARAM_HISTORY
 HISTORY_OUT_BASENAME      = "history",
 HISTORY_DEFAULT_TINTERVAL = 0.6D0,
 HISTORY_DEFAULT_TUNIT     = "SEC",
 HISTORY_DEFAULT_AVERAGE   = .false.,
 HISTORY_DATATYPE          = "REAL4",
/

#&HISTITEM item='DENS' /
#&HISTITEM item='MOMX' /
#&HISTITEM item='MOMY' /
#&HISTITEM item='MOMZ' /
#&HISTITEM item='RHOT' /
#&HISTITEM item='PRES' /
#&HISTITEM item='T'    /
#&HISTITEM item='U'    /
#&HISTITEM item='V'    /
#&HISTITEM item='W'    /
#&HISTITEM item='PT'   /

#&HISTITEM item='QV'   /
#&HISTITEM item='QC'   /
#&HISTITEM item='QR'   /
#&HISTITEM item='QI'   /
#&HISTITEM item='QS'   /
#&HISTITEM item='QG'   /
#&HISTITEM item='NC'   /
#&HISTITEM item='NR'   /
#&HISTITEM item='NI'   /
#&HISTITEM item='NS'   /
#&HISTITEM item='NG'   /

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
fpcoll -Ihwm,cpu -l0 -o Basic_Profile.txt -m 200000                        mpiexec $LPG $BIN/$EXE $EXE.cnf
#fpcoll -C -d pa1 -Usection=range,local_event_number=0,29,29,29,30,5,9,6    mpiexec $LPG $BIN/$EXE $EXE.cnf
#fpcoll -C -d pa2 -Usection=range,local_event_number=30,30,30,8,29,30,31,0  mpiexec $LPG $BIN/$EXE $EXE.cnf
#fpcoll -C -d pa3 -Usection=range,local_event_number=31,10,11,30,31,0,30,30 mpiexec $LPG $BIN/$EXE $EXE.cnf
#fpcoll -C -d pa4 -Usection=range,local_event_number=0,12,48,48,2,32,48,48  mpiexec $LPG $BIN/$EXE $EXE.cnf
#fpcoll -C -d pa5 -Usection=range,local_event_number=7,7,7,32,0,13,13,22    mpiexec $LPG $BIN/$EXE $EXE.cnf
#fpcoll -C -d pa6 -Usection=range,local_event_number=0,13,13,13,13,35,35,33 mpiexec $LPG $BIN/$EXE $EXE.cnf
#fpcoll -C -d pa7 -Usection=range,local_event_number=35,35,26,0,32,7,7,31   mpiexec $LPG $BIN/$EXE $EXE.cnf
echo "job ${RUNNAME} end     at " `date`

exit
