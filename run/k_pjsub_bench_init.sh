#! /bin/bash -x
#
# for K Computer
#
#PJM --rsc-list "node=2"
#PJM --rsc-list "elapse=01:00:00"
#PJM --rsc-list "node-mem=12Gi"
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=$PARALLEL
export LPG="lpgparm -s 4MB -d 4MB -h 256MB -t 4MB -p 4MB"

export HMDIR=/work/user0171/scale3

export BIN=/work/user0171/scale3/bin/K
export EXE=init_coldbubble

export OUTDIR=${HMDIR}/output/init_coldbubble

# Run Command
export RUN="fpcoll -Ibalance,call,cpu,hwm, -l20 -i20 -o Basic_Profile.txt -m 200000 mpiexec $LPG $BIN/$EXE init_coldbubble.cnf"

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
 PRC_NUM_X       = 2,
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
 GRID_OUT_BASENAME = 'grid_500m_70x100x220',
 GRID_DXYZ         = 40.D0,
 GRID_KMAX         = 220,
 GRID_IMAX         = 70,
 GRID_JMAX         = 100,
 GRID_BUFFER_DZ    =  4.D3,
 GRID_BUFFER_DX    =  0.D0,
 GRID_BUFFER_DY    =  0.D0,
 GRID_BUFFFACT     = 1.0D0,
/

&PARAM_ATMOS
 ATMOS_TYPE_DYN    = 'fent_fct',
 ATMOS_TYPE_PHY_TB = 'smagorinsky',
 ATMOS_TYPE_PHY_MP = 'NDW6',
 ATMOS_TYPE_PHY_RD = 'mstrnX',
/

&PARAM_ATMOS_VARS
 ATMOS_QTRC_NMAX            = 11,
 ATMOS_RESTART_OUT_BASENAME = 'init_coldbubble',
 ATMOS_RESTART_OUTPUT       = .true.,
/

&PARAM_MKEXP_COLDBUBBLE
 XC_BBL =   1.4D3,
 YC_BBL =   2.0D3,
 ZC_BBL =   4.0D3,
 XR_BBL =   5.0D2,
 YR_BBL =   5.0D2,
 ZR_BBL =   5.0D2,
/

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
$RUN > $OUTDIR/LOG 2>&1
echo "job ${RUNNAME} end     at " `date`

exit
