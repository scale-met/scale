#! /bin/bash -x

#BSUB -q as
#BSUB -n 6
#BSUB -J commtest
#BSUB -o commtest_log
#BSUB -e commtest_error

export HMDIR=${HOME}/scale3
export BIN=${HOME}/scale3/bin/Altix
export EXE=commtest

export OUTDIR=${HMDIR}/${EXE}

# Run Command
export MPIRUN="mpijob"

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
/

&PARAM_GRID
 GRID_OUT_BASENAME = 'grid_500m_120x120x25',
 GRID_DXYZ         = 10.D0,
 GRID_KMAX         = 25,
 GRID_IMAX         = 10,
 GRID_JMAX         = 10,
 GRID_BUFFER_DZ    =  0.D0,
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
 XC_BBL = 150.0D3,
 YC_BBL =  30.0D3,
 ZC_BBL =   4.0D3,
 XR_BBL =   3.0D3,
 YR_BBL =   3.0D3,
 ZR_BBL =   3.0D3,
/

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
$MPIRUN $BIN/$EXE ${EXE}.cnf > STDOUT 2>&1
echo "job ${RUNNAME} end     at " `date`

exit
