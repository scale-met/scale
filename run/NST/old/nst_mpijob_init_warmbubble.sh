#! /bin/bash -x
#
# for 3 team (NST) cluster
#
#BSUB -q as
#BSUB -n 12
#BSUB -J scale3
#BSUB -o scale3_log
#BSUB -e scale3_error


export OMP_NUM_THREADS=1

export HMDIR=/home/yashiro/scale3
export BIN=${HMDIR}/bin/NST-ifort
export EXE=init_warmbubble

export OUTDIR=${HMDIR}/output/init_warmbubble

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
 PRC_NUM_X       = 4,
 PRC_NUM_Y       = 3,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_TIME
 TIME_STARTDATE             = 2000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
/












&PARAM_GRID
 GRID_OUT_BASENAME = "grid_1000m_20x100x100",
 GRID_DXYZ         = 1000.D0,
 GRID_KMAX         = 20,
 GRID_IMAX         = 100,
 GRID_JMAX         = 100,
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

&PARAM_MKEXP_WARMBUBBLE
 EXT_TBBL =   3.D0,
 ZC_BBL   =   3.D3,
 XC_BBL   = 100.D3,
 YC_BBL   = 100.D3,
 ZR_BBL   =   3.D3,
 XR_BBL   =   6.D3,
 YR_BBL   =   6.D3,
/

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
mpijob $BIN/$EXE ${EXE}.cnf
echo "job ${RUNNAME} end     at " `date`

exit
