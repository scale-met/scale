#! /bin/bash -x
#
# for K Computer
#
#PJM --rsc-list "node=8x8"
#PJM --rsc-list "elapse=01:00:00"
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
export EXE=MPTB_NDW6

export OUTDIR=${HMDIR}/output/supercell-3

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
 PRC_NUM_X       = 8,
 PRC_NUM_Y       = 8,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_COMM
 COMM_dototalval = .true.,
/

&PARAM_TIME
 TIME_STARTDATE             = 2000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
 TIME_DURATION              = 1800.D0,
 TIME_DURATION_UNIT         = 'SEC',
 TIME_DT                    = 1.2D0,
 TIME_DT_UNIT               = 'SEC',
 TIME_DT_ATMOS_DYN          = 0.6D0,
 TIME_DT_ATMOS_DYN_UNIT     = 'SEC',
 TIME_NSTEP_ATMOS_DYN       = 2,
 TIME_DT_ATMOS_PHY_MP       = 1.2D0,
 TIME_DT_ATMOS_PHY_MP_UNIT  = 'SEC',
 TIME_DT_ATMOS_PHY_TB       = 1.2D0,
 TIME_DT_ATMOS_PHY_TB_UNIT  = "SEC",
 TIME_DT_ATMOS_RESTART      = 900.D0,
 TIME_DT_ATMOS_RESTART_UNIT = 'SEC',
/

&PARAM_GRID
 GRID_OUT_BASENAME = "",
 GRID_DXYZ         = 400.D0,
 GRID_KMAX         = 50,
 GRID_IMAX         = 50,
 GRID_JMAX         = 50,
 GRID_BUFFER_DZ    = 8.0D3,
 GRID_BUFFFACT     = 1.0D0,
/

&PARAM_ATMOS
 ATMOS_TYPE_DYN    = "fent_fct",
 ATMOS_TYPE_PHY_TB = "smagorinsky",
 ATMOS_TYPE_PHY_MP = "NDW6",
 ATMOS_TYPE_PHY_RD = "mstrnX",
/

&PARAM_ATMOS_VARS
 ATMOS_QTRC_NMAX              = 11,
 ATMOS_RESTART_IN_BASENAME    = '${HMDIR}/output/supercell-2/init_supercell_63072003600.000',
 ATMOS_RESTART_OUTPUT         = .true.,
 ATMOS_RESTART_OUT_BASENAME   = 'init_supercell',
/

&PARAM_ATMOS_REFSTATE
 ATMOS_REFSTATE_TEMP_SFC    =   300.D0     
/

&PARAM_ATMOS_BOUNDARY
 ATMOS_BOUNDARY_OUT_BASENAME = "boundary",
 ATMOS_BOUNDARY_VALUE_VELX   =  0.D0,
 ATMOS_BOUNDARY_TAUZ         = 10.D0,
/

&PARAM_ATMOS_DYN
 ATMOS_DYN_NUMERICAL_DIFF = 1.D-2,
/

#&PARAM_ATMOS_PHY_MP
# doreport_tendency = .true.,
#/

&NM_MP_NDW6_INIT
! OPT_DEBUG     = .true.,
! OPT_DEBUG_TEM = .true.,
/
&PARAM_OCEAN
 OCEAN_TYPE = 'FIXEDSST',
/

&PARAM_HISTORY
 HISTORY_OUT_BASENAME      = 'history',
 HISTORY_DEFAULT_TINTERVAL = 60.D0,
 HISTORY_DEFAULT_TUNIT     = 'SEC',
 HISTORY_DEFAULT_AVERAGE   = .false.,
 HISTORY_DATATYPE          = 'REAL4',
/

&HISTITEM item='DENS' /
#&HISTITEM item='MOMX' /
#&HISTITEM item='MOMY' /
#&HISTITEM item='MOMZ' /
#&HISTITEM item='RHOT' /
#&HISTITEM item='PRES' /
&HISTITEM item='T'    /
&HISTITEM item='U'    /
&HISTITEM item='V'    /
&HISTITEM item='W'    /
&HISTITEM item='PT'   /

&HISTITEM item='QTOT' /
&HISTITEM item='QV'   /
&HISTITEM item='QC'   /
&HISTITEM item='QR'   /
&HISTITEM item='QI'   /
&HISTITEM item='QS'   /
&HISTITEM item='QG'   /
&HISTITEM item='NC'   /
&HISTITEM item='NR'   /
&HISTITEM item='NI'   /
&HISTITEM item='NS'   /
&HISTITEM item='NG'   /

#&HISTITEM item='dqprcp' /
#&HISTITEM item='dqcond' /
#&HISTITEM item='dqevap' /
#&HISTITEM item='dqauto' /
#&HISTITEM item='dqcoll' /

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
fpcoll -Ihwm,cpu -l0 -o Basic_Profile.txt -m 200000                        mpiexec $LPG $BIN/$EXE $EXE.cnf
echo "job ${RUNNAME} end     at " `date`

exit
