#! /bin/bash -x

export HMDIR=~/GCMresults/sol/latest/output
export BIN=~/Dropbox/Inbox/scale3/bin/MacOSX-ifort
export EXE=SmagorinskyAndSurface

export OUTDIR=${HMDIR}/data/SmagorinskyAndSurface

# Run Command
export MPIRUN="/usr/local/mpich213/bin/mpiexec -np 6 -f /Users/yashiro/libs/mpilib/machines_local"

mkdir -p $OUTDIR
cd $OUTDIR

########################################################################
cat << End_of_SYSIN > ${OUTDIR}/${EXE}.cnf

#####
#
# Scale3 benchmark configuration
#
#####

&PARAM_PRC
 PRC_NUM_X       = 3,
 PRC_NUM_Y       = 2,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_CONST
/

&PARAM_TIME
 TIME_STARTDATE             = 2000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
 TIME_DURATION              = 40.D0,
 TIME_DURATION_UNIT         = 'SEC',
 TIME_DT                    = 0.4D0,
 TIME_DT_UNIT               = 'SEC',
 TIME_DT_ATMOS_DYN          = 0.4D0,
 TIME_DT_ATMOS_DYN_UNIT     = 'SEC',
 TIME_NSTEP_ATMOS_DYN       = 1,
 TIME_DT_ATMOS_PHY_TB       = 0.4D0,
 TIME_DT_ATMOS_PHY_TB_UNIT  = 'SEC',
 TIME_DT_ATMOS_RESTART      = 40.D0,
 TIME_DT_ATMOS_RESTART_UNIT = 'SEC',
/

&PARAM_GRID
 GRID_OUT_BASENAME = 'grid_500m_120x120x25',
 GRID_DXYZ         = 500.D0,
 GRID_KMAX         = 25,
 GRID_IMAX         = 120,
 GRID_JMAX         = 120,
 GRID_BUFFER_DZ    = 4.0D3,
 GRID_BUFFER_DX    = 0.0D0,
 GRID_BUFFER_DY    = 0.0D0,
 GRID_BUFFFACT     = 1.0D0,
/

&PARAM_ATMOS
 ATMOS_TYPE_DYN    = 'fent_fct',
 ATMOS_TYPE_PHY_TB = 'smagorinsky',
 ATMOS_TYPE_PHY_MP = 'NDW6',
 ATMOS_TYPE_PHY_RD = 'mstrnX',
/

&PARAM_ATMOS_VARS
 ATMOS_QTRC_NMAX              = 11,
 ATMOS_RESTART_IN_BASENAME    = '${HMDIR}/data/init_coldbubble_lores/init_coldbubble_63072000000.000',
 ATMOS_RESTART_OUTPUT         = .true.,
 ATMOS_RESTART_OUT_BASENAME   = 'check_coldbubble',
 ATMOS_RESTART_CHECK          = .false.,
 ATMOS_RESTART_CHECK_BASENAME = '${HMDIR}/data/coldbubble_lores/check_coldbubble_63072000040.000',
/

&PARAM_ATMOS_REFSTATE
 ATMOS_REFSTATE_TEMP_SFC    =   300.D0     
/
# ATMOS_REFSTATE_OUT_BASENAME = 'refstate',

&PARAM_ATMOS_BOUNDARY
 ATMOS_BOUNDARY_TAUZ = 75.D0,
/

&PARAM_ATMOS_DYN
 ATMOS_DYN_NUMERICAL_DIFF = 1.D-3,
/

&PARAM_OCEAN
 OCEAN_TYPE = 'FIXEDSST',
/

&PARAM_OCEAN_VARS
 OCEAN_SST = 293.15D0,
/

&PARAM_HISTORY
 HISTORY_OUT_BASENAME      = 'history',
 HISTORY_DEFAULT_TINTERVAL = 4.D0,
 HISTORY_DEFAULT_TUNIT     = 'SEC',
 HISTORY_DEFAULT_AVERAGE   = .false.,
 HISTORY_DATATYPE          = 'REAL4',
/

&HISTITEM item='DENS' /
&HISTITEM item='MOMX' /
&HISTITEM item='MOMY' /
&HISTITEM item='MOMZ' /
&HISTITEM item='RHOT' /
&HISTITEM item='MOMX_t_tb' /
&HISTITEM item='MOMY_t_tb' /
&HISTITEM item='MOMZ_t_tb' /
&HISTITEM item='RHOT_t_tb' /
&HISTITEM item='PRES' /
&HISTITEM item='U'    /
&HISTITEM item='V'    /
&HISTITEM item='W'    /
&HISTITEM item='T'    /
&HISTITEM item='PT'   /

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
$MPIRUN $BIN/$EXE ${EXE}.cnf > STDOUT 2>&1
echo "job ${RUNNAME} end     at " `date`

exit
