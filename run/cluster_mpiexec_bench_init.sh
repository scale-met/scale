#! /bin/bash -x
#
# for AICS.common
#
#$ -pe  mpi 1
#$ -q   large
#$ -cwd

export RLIMIT_MEMLOCK=65536
ulimit -l 65536
. /opt/intel/ics/2011.05.23/ictvars.sh

export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_HYDRA_BOOTSTRAP_EXEC=/usr/bin/rsh

export HMDIR=/home/yashiro/scale3
export BIN=${HMDIR}/bin/${SCALE_SYS}
export EXE=init_coldbubble

export OUTDIR=${HMDIR}/output/init_coldbubble

mkdir -p ${OUTDIR}
cd ${OUTDIR}

########################################################################
cat << End_of_SYSIN > ${OUTDIR}/${EXE}.cnf

#####
#
# Scale3 init_coldbubble configulation
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
 GRID_OUT_BASENAME = "grid_20m_336x63x63",
 GRID_DXYZ         = 20.D0,
 GRID_KMAX         = 336,
 GRID_IMAX         = 63,
 GRID_JMAX         = 63,
 GRID_BUFFER_DZ    = 6.0D3,
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
 ATMOS_RESTART_OUT_BASENAME   = "init_coldbubble",
/

&PARAM_MKEXP_COLDBUBBLE
 ZC_BBL =  6.0D2,
 XC_BBL =  6.0D2,
 YC_BBL =  6.0D2,
 ZR_BBL =  2.0D2,
 XR_BBL =  2.0D2,
 YR_BBL =  2.0D2,
/

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
mpiexec.hydra -hostfile $TMPDIR/machines -n 1 -perhost 1 $BIN/$EXE ${EXE}.cnf
echo "job ${RUNNAME} end     at " `date`

exit
