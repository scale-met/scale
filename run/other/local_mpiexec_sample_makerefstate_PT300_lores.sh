#! /bin/bash -x

export HMDIR=~/GCMresults/sol/latest/output
export BIN=~/Dropbox/Inbox/scale3/bin/MacOSX-ifort
export EXE=makerefstate

export OUTDIR=${HMDIR}/data/makerefstate

# Run Command
export MPIRUN="/usr/local/mpich213/bin/mpiexec -np 1 -f /Users/yashiro/libs/mpilib/machines_local"

mkdir -p $OUTDIR
cd $OUTDIR

########################################################################
cat << End_of_SYSIN > ${OUTDIR}/${EXE}.cnf

#####
#
# Scale3 init_coldbubble configulation
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
/

&PARAM_GRID
 GRID_DXYZ         = 500.D0,
 GRID_KMAX         = 25,
 GRID_IMAX         = 120,
 GRID_JMAX         = 120,
 GRID_BUFFER_DZ    = 4.0D3,
 GRID_BUFFER_DX    = 0.0D0,
 GRID_BUFFER_DY    = 0.0D0,
 GRID_BUFFFACT     = 1.0D0,
/

&PARAM_ATMOS_REFSTATE
 ATMOS_REFSTATE_OUT_BASENAME = 'refstate_PT300_500m_z25',
 ATMOS_REFSTATE_TYPE         = 'UNIFORM',
 ATMOS_REFSTATE_POTT_UNIFORM = 300.D0
/

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
$MPIRUN $BIN/$EXE ${EXE}.cnf > STDOUT 2>&1
echo "job ${RUNNAME} end     at " `date`

exit
