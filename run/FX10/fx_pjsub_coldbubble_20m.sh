#! /bin/bash -x
#
# for OAKLEAF-FX
#
#PJM -L "rscgrp=short"
#PJM -L "node=1x1"
#PJM -L "elapse=00:10:00"
#PJM -j

export PARALLEL=8
export OMP_NUM_THREADS=$PARALLEL
export fu08bf=10

export HMDIR=/home/c18001/scale3
export BIN=${HMDIR}/bin/FX10
export EXE=scale3_336x63x63_ndw6_

export OUTDIR=${HMDIR}/output/coldbubble_20m

mkdir -p ${OUTDIR}
cd ${OUTDIR}
rm -rf ./prof

########################################################################
cat << End_of_SYSIN > ${OUTDIR}/${EXE}.cnf

#####
#
# Scale3 benchmark configuration
#
#####

&PARAM_IO
 IO_LOG_SUPPRESS = .false.,
 IO_LOG_ALLNODE  = .false.,
/

&PARAM_PRC
 PRC_NUM_X       = 1,
 PRC_NUM_Y       = 1,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_CONST
/

&PARAM_TIME
 TIME_STARTDATE             = 2000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
 TIME_DURATION              = 30.D0,
 TIME_DURATION_UNIT         = "SEC",
 TIME_DT                    = 0.6D0,
 TIME_DT_UNIT               = "SEC",
 TIME_DT_ATMOS_DYN          = 0.03D0,
 TIME_DT_ATMOS_DYN_UNIT     = "SEC",
 TIME_NSTEP_ATMOS_DYN       = 20,
 TIME_DT_ATMOS_PHY_TB       = 0.6D0,
 TIME_DT_ATMOS_PHY_TB_UNIT  = "SEC",
 TIME_DT_ATMOS_PHY_MP       = 0.6D0,
 TIME_DT_ATMOS_PHY_MP_UNIT  = "SEC",
 TIME_DT_ATMOS_PHY_RD       = 0.6D0,
 TIME_DT_ATMOS_PHY_RD_UNIT  = "SEC",
 TIME_DT_ATMOS_RESTART      = 15.D0,
 TIME_DT_ATMOS_RESTART_UNIT = "SEC",
 TIME_DT_OCEAN              = 15.D0,
 TIME_DT_OCEAN_UNIT         = "MIN",
/

&PARAM_GRID
 GRID_IN_BASENAME  = "",
 GRID_OUT_BASENAME = "",
/

&PARAM_GEOMETRICS
 GEOMETRICS_startlonlat  = 120.D0, 30.D0,
 GEOMETRICS_rotation     = 10.D0,
 GEOMETRICS_OUT_BASENAME = "",
/

&PARAM_COMM
 COMM_vsize_max       = 40,
 COMM_total_doreport  = .true.,
 COMM_total_globalsum = .true.,
/

&PARAM_ATMOS
 ATMOS_TYPE_DYN    = "fent_pdfct",
 ATMOS_TYPE_PHY_TB = "smagorinsky",
 ATMOS_TYPE_PHY_MP = "NDW6",
 ATMOS_TYPE_PHY_RD = "mstrnX",
/

&PARAM_ATMOS_VARS
 ATMOS_RESTART_IN_BASENAME      = "${HMDIR}/output/init_coldbubble_20m/init_coldbubble_63072000000.000",
 ATMOS_RESTART_IN_ALLOWMISSINGQ = .false.,
 ATMOS_RESTART_OUTPUT           = .true.,
 ATMOS_RESTART_OUT_BASENAME     = "init_coldbubble",
 ATMOS_RESTART_CHECK            = .false.,
 ATMOS_RESTART_CHECK_BASENAME   = "restart_check",
 ATMOS_RESTART_CHECK_CRITERION  =  1.D-6
/

&PARAM_ATMOS_REFSTATE
 ATMOS_REFSTATE_IN_BASENAME  = "",
 ATMOS_REFSTATE_OUT_BASENAME = "",
 ATMOS_REFSTATE_TYPE         = "UNIFORM",
 ATMOS_REFSTATE_POTT_UNIFORM = 300.D0
 ATMOS_REFSTATE_TEMP_SFC     = 300.D0     
/

&PARAM_ATMOS_BOUNDARY
 ATMOS_BOUNDARY_OUT_BASENAME = "boundary",
 ATMOS_BOUNDARY_VALUE_VELX   =  0.D0,
 ATMOS_BOUNDARY_TAUZ         = 10.D0,
/

&PARAM_ATMOS_DYN
 ATMOS_DYN_NUMERICAL_DIFF = 1.D-2,
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

&HISTITEM item='U'    /
&HISTITEM item='V'    /
&HISTITEM item='W'    /
&HISTITEM item='PT'   /
#&HISTITEM item='PRES' /
#&HISTITEM item='T'    /

&HISTITEM item='VOR'  /
&HISTITEM item='ENGP' /
&HISTITEM item='ENGK' /
&HISTITEM item='ENGI' /

&HISTITEM item='QV'   /
#&HISTITEM item='QTOT' /
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

&HISTITEM item='TKE'  /
&HISTITEM item='NU'   /
#&HISTITEM item='Pr'   /
#&HISTITEM item='Ri'   /

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
fipp -C -Ihwm -d prof mpiexec $BIN/$EXE $EXE.cnf
echo "job ${RUNNAME} end     at " `date`

exit
