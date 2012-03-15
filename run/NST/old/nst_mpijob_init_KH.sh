#! /bin/bash -x
#
# for 3 team (NST) cluster
#
#BSUB -q cl
#BSUB -n 96
#BSUB -J scale3
#BSUB -o scale3_log
#BSUB -e scale3_error

export OMP_NUM_THREADS=1

export HMDIR=/home/yashiro/scale3
export BIN=${HMDIR}/bin/NST-ifort
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
 PRC_NUM_X       = 16,
 PRC_NUM_Y       = 6,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_TIME
 TIME_STARTDATE             = 2000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
/












&PARAM_GRID
 GRID_OUT_BASENAME = "grid_40m_220x25x25",
 GRID_DXYZ         = 40.D0,
 GRID_KMAX         = 220,
 GRID_IMAX         = 25,
 GRID_JMAX         = 25,
 GRID_BUFFER_DZ    = 4.0D3,
 GRID_BUFFFACT     = 1.1D0,
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
mpijob $BIN/$EXE ${EXE}.cnf > STDOUT 2>&1
echo "job ${RUNNAME} end     at " `date`

exit
