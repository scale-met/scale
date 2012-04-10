#! /bin/bash -x
#
# for Local machine (MacOSX 2core+gfortran+MPICH2)
#
export HMDIR=~/Desktop/scale3
export BIN=~/Dropbox/Inbox/scale3/bin/MacOSX-gfortran
export EXE=init_supercell

export OUTDIR=${HMDIR}/output/init_supercell

mkdir -p ${OUTDIR}
cd ${OUTDIR}

########################################################################
cat << End_of_SYSIN > ${OUTDIR}/${EXE}.cnf

#####
#
# Scale3 init_supercell configulation
#
#####

&PARAM_PRC
 PRC_NUM_X       = 1,
 PRC_NUM_Y       = 1,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_TIME
 TIME_STARTDATE             = 2000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
/












&PARAM_GRID
 GRID_OUT_BASENAME = "grid_400m_46x168x168",
 GRID_DXYZ         = 400.D0,
 GRID_KMAX         = 46,
 GRID_IMAX         = 168,
 GRID_JMAX         = 168,
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
 ATMOS_QTRC_NMAX            = 11,
 ATMOS_RESTART_OUTPUT       = .true.,
 ATMOS_RESTART_OUT_BASENAME = 'init_supercell',
/

&PARAM_MKEXP_SUPERCELL
 ENV_IN_SOUNDING_file = '${HMDIR}/data/supercell/input_sounding.txt',
 ZC_BBL =  1.5D3,
 XC_BBL = 67.2D2,
 YC_BBL = 67.2D2,
 ZR_BBL =  1.0D3,
 XR_BBL =  8.0D3,
 YR_BBL =  8.0D3,
/

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
/opt/local/bin/mpiexec -np 1 -f /Users/yashiro/machines_local $BIN/$EXE ${EXE}.cnf > STDOUT 2>&1
echo "job ${RUNNAME} end     at " `date`

exit
