#! /bin/bash -x
#
# for K Computer
#
#PJM --rsc-list "node=24x4"
#PJM --rsc-list "elapse=00:10:00"
#PJM --rsc-list "node-mem=10Gi"
#PJM -s
#
. /work/system/Env_base_1.2.0-02
#
export PARALLEL=8
export OMP_NUM_THREADS=$PARALLEL
export LPG="lpgparm -s 32MB -d 32MB -h 32MB -t 32MB -p 32MB"
export fu08bf=1

export HMDIR=/work/scratch/user0171/scale3
export BIN=${HMDIR}/bin/K
export EXE=init_turbdyn

export OUTDIR=${HMDIR}/output/init_KH

mkdir -p ${OUTDIR}
cd ${OUTDIR}

########################################################################
cat << End_of_SYSIN > ${OUTDIR}/${EXE}.cnf

#####
#
# Scale3 init_warmbubble configulation
#
#####

&PARAM_PRC
 PRC_NUM_X       = 24,
 PRC_NUM_Y       = 4,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_TIME
 TIME_STARTDATE             = 2000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
/












&PARAM_GRID
 GRID_OUT_BASENAME = "grid_20m_200x25x25",
 GRID_DXYZ         = 20.D0,
 GRID_KMAX         = 200,
 GRID_IMAX         = 25,
 GRID_JMAX         = 25,
 GRID_BUFFER_DZ    = 0.0D3,
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
 ATMOS_RESTART_OUTPUT         = .true.,
 ATMOS_RESTART_OUT_BASENAME   = "init_KH",
/

&PARAM_MKEXP_TURBDYN
/

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
fpcoll -Ihwm,cpu -l0 -o Basic_Profile.txt -m 200000 mpiexec $LPG $BIN/$EXE $EXE.cnf
echo "job ${RUNNAME} end     at " `date`

exit
