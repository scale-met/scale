#!/bin/bash
#$ -pe mpi 1
#$ -q large
#$ -cwd

export HMDIR=/home/shimosaka/CLIMATE/CMP3/scale3
export BIN=${HMDIR}/bin/${SCALE_SYS}
export EXE=scale3

export OUTDIR=${HMDIR}/data/coldbubble

# Run Command
export RLIMIT_MEMLOCK=65536
ulimit -l 65536
. /opt/intel/ics/2011.05.23/ictvars.sh
export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_HYDRA_BOOTSTRAP_EXEC=/usr/bin/rsh

export MPIRUN="mpiexec.hydra -hostfile $TMPDIR/machines -n 1 -perhost 1"
#export MPIRUN="mpiexec.hydra -n 1"

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
 TIME_DURATION              = 6.0D0,
 TIME_DURATION_UNIT         = 'SEC',
 TIME_DT                    = 0.6D0,
 TIME_DT_UNIT               = 'SEC',
 TIME_DT_ATMOS_DYN          = 0.06D0,
 TIME_DT_ATMOS_DYN_UNIT     = 'SEC',
 TIME_NSTEP_ATMOS_DYN       = 10,
 TIME_DT_ATMOS_RESTART      = 6.D0,
 TIME_DT_ATMOS_RESTART_UNIT = 'SEC',
/

&PARAM_GRID
 GRID_OUT_BASENAME = 'grid_40m_70x100x220',
 GRID_DXYZ         = 40.D0,
 GRID_KMAX         = 220,
 GRID_IMAX         = 70,
 GRID_JMAX         = 100,
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
 ATMOS_RESTART_IN_BASENAME    = '${HMDIR}/data/init_coldbubble/init_coldbubble_63072000000.000',
 ATMOS_RESTART_OUTPUT         = .false.,
 ATMOS_RESTART_OUT_BASENAME   = 'check_coldbubble',
 ATMOS_RESTART_CHECK          = .true.,
 ATMOS_RESTART_CHECK_BASENAME = '${HMDIR}/data/coldbubble/check_coldbubble_63072000006.000',
/

&PARAM_ATMOS_REFSTATE
 ATMOS_REFSTATE_OUT_BASENAME = 'refstate',
 ATMOS_REFSTATE_TEMP_SFC    =   300.D0     
/

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
 HISTORY_DEFAULT_TINTERVAL = 0.6D0,
 HISTORY_DEFAULT_TUNIT     = 'MSEC',
 HISTORY_DEFAULT_AVERAGE   = .false.,
 HISTORY_DATATYPE          = 'REAL4',
/

#&HISTITEM item='DENS' /
#&HISTITEM item='MOMX' /
#&HISTITEM item='MOMY' /
#&HISTITEM item='MOMZ' /
#&HISTITEM item='RHOT' /
#&HISTITEM item='PRES' /
#&HISTITEM item='U'    /
#&HISTITEM item='V'    /
#&HISTITEM item='W'    /
#&HISTITEM item='T'    /
#&HISTITEM item='PT'   /

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
${MPIRUN} ${BIN}/${EXE} ${OUTDIR}/${EXE}.cnf > STDOUT 2>&1
echo "job ${RUNNAME} end     at " `date`

exit
