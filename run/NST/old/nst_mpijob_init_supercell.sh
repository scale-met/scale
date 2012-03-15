#! /bin/bash -x
#
# for 3 team (NST) cluster
#
#BSUB -q cl
#BSUB -n 64
#BSUB -J scale3
#BSUB -o scale3_log
#BSUB -e scale3_error
#BSUB -R "span[ptile=8]"

export OMP_NUM_THREADS=1

export HMDIR=/home/yashiro/scale3
export BIN=${HMDIR}/bin/NST-ifort
export EXE=init_supercell

export OUTDIR=${HMDIR}/output/init_supercell3

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
 PRC_NUM_X       = 8,
 PRC_NUM_Y       = 8,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_TIME
 TIME_STARTDATE             = 2000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
/












&PARAM_GRID
 GRID_OUT_BASENAME = 'grid_400m_50x50x50',
 GRID_DXYZ         = 400.D0,
 GRID_KMAX         = 50,
 GRID_IMAX         = 50,
 GRID_JMAX         = 50,
 GRID_BUFFER_DZ    = 8.0D3,
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
 ENV_IN_SOUNDING_file = '${HMDIR}/data/supercell/input_sounding.txt',
 EXT_TBBL = 3.D0
 ZC_BBL =  1.5D3,
 XC_BBL = 80.0D3,
 YC_BBL = 80.0D3,
 ZR_BBL =  3.0D3,
 XR_BBL = 20.0D3,
 YR_BBL = 20.0D3,
/

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
mpijob $BIN/$EXE ${EXE}.cnf
echo "job ${RUNNAME} end     at " `date`

exit
