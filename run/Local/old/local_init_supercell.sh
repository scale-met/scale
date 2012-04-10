#! /bin/bash -x
#
# for Local machine (MacOSX 8core+ifort+MPICH2)
#
export HMDIR=~/GCMresults/sol/latest
export BIN=~/Dropbox/Inbox/scale3/bin/${SCALE_SYS}
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
 PRC_NUM_X       = 2,
 PRC_NUM_Y       = 2,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_TIME
 TIME_STARTDATE             = 2000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
/












&PARAM_GRID
 GRID_OUT_BASENAME = 'grid_400m_50x100x100',
 GRID_DXYZ         = 400.D0,
 GRID_KMAX         = 50,
 GRID_IMAX         = 100,
 GRID_JMAX         = 100,
 GRID_BUFFER_DZ    = 4.0D3,
 GRID_BUFFER_DX    = 0.0D0,
 GRID_BUFFER_DY    = 0.0D0,
 GRID_BUFFFACT     = 1.0D0,
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
 ENV_IN_SOUNDING_file = '~/Dropbox/Inbox/scale3/data/supercell/input_sounding.txt',
 ZC_BBL =  1.5D3,
 XC_BBL = 40.0D3,
 YC_BBL = 40.0D3,
 ZR_BBL =  1.0D3,
 XR_BBL =  4.0D3,
 YR_BBL =  4.0D3,
/

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
/usr/local/mpich213/bin/mpiexec -np 4 -f /Users/yashiro/libs/mpilib/machines_local $BIN/$EXE ${EXE}.cnf > STDOUT 2>&1
echo "job ${RUNNAME} end     at " `date`

exit
