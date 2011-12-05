#! /bin/bash -x
#
# for K Computer
#
#PJM --rsc-list "node=1024"
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
export EXE=init_coldbubble

export OUTDIR=${HMDIR}/output/coldbubble_32x32

# Run Command
export RUN="fpcoll -Ihwm -Srange -o Basic_Profile.txt -m 200000 mpiexec $LPG $BIN/$EXE init_coldbubble.cnf"

mkdir -p $OUTDIR
cd $OUTDIR

########################################################################
cat << End_of_SYSIN > ${OUTDIR}/init_coldbubble.cnf

#####
#
# Scale3 init_coldbubble configulation
#
#####

&PARAM_TIME
 TIME_STARTDATE             =  2000, 1, 1, 0, 0, 0,
 TIME_STARTMS               =  0.D0,
/









&PARAM_PRC
 PRC_NUM_X = 32,
 PRC_NUM_Y = 32,
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
 ATMOS_RESTART_OUT_BASENAME = 'init_coldbubble',
/

&PARAM_MKEXP_COLDBUBBLE
 XC_BBL = 0.8D3,
 YC_BBL = 0.8D3,
 ZC_BBL = 4.0D3,
 XR_BBL = 0.5D3,
 YR_BBL = 0.5D3,
 ZR_BBL = 1.5D3,
/

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
$RUN > $OUTDIR/LOG 2>&1
echo "job ${RUNNAME} end     at " `date`

exit
