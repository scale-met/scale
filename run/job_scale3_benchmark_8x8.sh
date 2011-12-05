#! /bin/bash -x
#
# for K Computer
#
#PJM --rsc-list "node=64"
#PJM --rsc-list "elapse=01:00:00"
#PJM --rsc-list "node-mem=12Gi"
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=$PARALLEL
export LPG="lpgparm -s 4MB -d 4MB -h 256MB -t 4MB -p 4MB"

export HMDIR=/work/user0171

export BIN=/work/user0171/scale3/bin
export EXE=scale3

export OUTDIR=${HMDIR}/output/scale3_8x8

# Run Command
export RUN="fpcoll -Ihwm -o Basic_Profile.txt -i100 -m 200000 mpiexec $LPG $BIN/$EXE scale3.cnf"

mkdir -p $OUTDIR
cd $OUTDIR

########################################################################
cat << End_of_SYSIN > ${OUTDIR}/scale3.cnf

#####
#
# Scale3 benchmark configuration
#
#####

&PARAM_TIME
 TIME_STARTDATE             =  2000, 1, 1, 0, 0, 0,
 TIME_STARTMS               =  0.D0,
 TIME_DURATION              =  30D0,
 TIME_DURATION_UNIT         = 'SEC',
 TIME_DT                    =  60.D0,
 TIME_DT_UNIT               = 'MSEC',
 TIME_DT_ATMOS_DYN          =  60.D0,
 TIME_DT_ATMOS_DYN_UNIT     = 'MSEC',
 TIME_DT_ATMOS_RESTART      =  60.D0,
 TIME_DT_ATMOS_RESTART_UNIT = 'SEC',
/

&PARAM_PRC
 PRC_NUM_X = 8,
 PRC_NUM_Y = 8,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_GRID
 GRID_DX =   40.D0,
 GRID_DY =   40.D0,
 GRID_DZ =   40.D0,
/

&PARAM_ATMOS
 ATMOS_TYPE_DYN    = 'fent_fct',
 ATMOS_TYPE_PHY_TB = 'smagorinsky',
 ATMOS_TYPE_PHY_MP = 'NDW6',
 ATMOS_TYPE_PHY_RD = 'mstrnX',
/

&PARAM_ATMOS_VARS
 ATMOS_QTRC_NMAX = 11,
 ATMOS_RESTART_IN_BASENAME = '../coldbubble_8x8/init_coldbubble_63072000000.000',
/

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
$RUN > $OUTDIR/LOG 2>&1
echo "job ${RUNNAME} end     at " `date`

exit
