#! /bin/bash -x
#
# for Local machine (MacOSX 8core+ifort+MPICH2)
#
export HMDIR=~/GCMresults/sol/latest
export BIN=~/Dropbox/Inbox/scale3/bin/${SCALE_SYS}
export EXE=init_warmbubble

export OUTDIR=${HMDIR}/output/init_warmbubble

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
 PRC_NUM_X       = 3,
 PRC_NUM_Y       = 2,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_TIME
 TIME_STARTDATE             = 2000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
/












&PARAM_GRID
 GRID_OUT_BASENAME = "grid_1000m_19x126x126",
 GRID_DXYZ         = 1000.D0,
 GRID_KMAX         = 19,
 GRID_IMAX         = 126,
 GRID_JMAX         = 126,
 GRID_BUFFER_DZ    = 5.0D3,
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
 ATMOS_RESTART_OUT_BASENAME   = "init_warmbubble",
/

&PARAM_ATMOS_REFSTATE
 ATMOS_REFSTATE_TEMP_SFC = 300.D0     
/

&PARAM_MKEXP_WARMBUBBLE
 ZC_BBL =    6.D3,
 XC_BBL =  126.D3,
 YC_BBL =  126.D3,
 ZR_BBL =    4.D3,
 XR_BBL =   20.D3,
 YR_BBL =   20.D3,
/

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
/usr/local/mpich213/bin/mpiexec -np 6 -f /Users/yashiro/libs/mpilib/machines_local $BIN/$EXE ${EXE}.cnf > STDOUT 2>&1
echo "job ${RUNNAME} end     at " `date`

exit
