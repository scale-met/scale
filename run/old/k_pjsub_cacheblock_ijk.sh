#! /bin/bash -x
#
# for K Computer
#
#PJM --rsc-list "node=1x1"
#PJM --rsc-list "elapse=01:00:00"
#PJM --rsc-list "node-mem=12Gi"
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=$PARALLEL
export LPG="lpgparm -s 32MB -d 32MB -h 32MB -t 32MB -p 32MB"

export HMDIR=/work/user0171/scale3
export BIN=/work/user0171/scale3/bin/K
export EXE=cacheblock

export OUTDIR=${HMDIR}/output/cacheblock_ijk_1x1

mkdir -p $OUTDIR
cd $OUTDIR

echo "job ${RUNNAME} started at " `date`

for i in 1 3 5 7 9
do
for j in 7 23 39
do
for k in 220 110 55
do

########################################################################
cat << End_of_SYSIN > ${OUTDIR}/${EXE}.cnf

&PARAM_PRC
 PRC_NUM_X       = 1,
 PRC_NUM_Y       = 1,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/
&PARAM_TIME
 TIME_STARTDATE             = 2000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
 TIME_DURATION              = 6.D0,
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
 GRID_OUT_BASENAME = '',
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
 ATMOS_RESTART_IN_BASENAME    = '${HMDIR}/output/init_warmbubble/init_warmbubble_63072000000.000',
 ATMOS_RESTART_OUTPUT         = .false.,
 ATMOS_RESTART_CHECK          = .false.,
/
&PARAM_ATMOS_REFSTATE ATMOS_REFSTATE_TEMP_SFC = 300.D0, /
&PARAM_OCEAN OCEAN_TYPE = 'FIXEDSST', /
&PARAM_OCEAN_VARS OCEAN_SST = 293.15D0, /

&PARAM_ATMOS_DYN
 ATMOS_DYN_NUMERICAL_DIFF = 1.D-3,
 KBLOCK = ${k},
 IBLOCK = ${i},
 JBLOCK = ${j},
/

&PARAM_HISTORY
 HISTORY_OUT_BASENAME      = 'history',
 HISTORY_DEFAULT_TINTERVAL = 0.6D0,
 HISTORY_DEFAULT_TUNIT     = 'SEC',
 HISTORY_DEFAULT_AVERAGE   = .false.,
 HISTORY_DATATYPE          = 'REAL4',
/

End_of_SYSIN
########################################################################

fpcoll -Ihwm,cpu -l0 -o Basic_Profile_k${k}xi${i}xj${j}.txt -m 200000 mpiexec $LPG $BIN/$EXE $EXE.cnf > $OUTDIR/LOG 2>&1

done
done
done

echo "job ${RUNNAME} finished at " `date`

exit
