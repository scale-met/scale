#! /bin/bash -x
#
# for AICS.common
#
#$ -pe  mpi 1
#$ -q   small
#$ -cwd

export RLIMIT_MEMLOCK=65536
ulimit -l 65536
. /opt/intel/ics/2011.05.23/ictvars.sh

export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_HYDRA_BOOTSTRAP_EXEC=/usr/bin/rsh

export HMDIR=/home/yashiro/scale3
export BIN=${HMDIR}/bin/${SCALE_SYS}
export EXE=scale3_init_336x63x63_ndw6_

export OUTDIR=${HMDIR}/output/init_coldbubble

mkdir -p ${OUTDIR}
cd ${OUTDIR}

########################################################################
cat << End_of_SYSIN > ${OUTDIR}/${EXE}.cnf

#####
#
# Scale3 mkinit configulation
#
#####

&PARAM_PRC
 PRC_NUM_X       = 1,
 PRC_NUM_Y       = 1,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_GRID
 GRID_OUT_BASENAME = "grid_020m_336x63x63",
/

&PARAM_GEOMETRICS
 GEOMETRICS_startlonlat  = 120.D0, 30.D0,
 GEOMETRICS_rotation     = 10.D0,
 GEOMETRICS_OUT_BASENAME = "",
/

&PARAM_COMM
 COMM_total_doreport  = .true.,
 COMM_total_globalsum = .true.,
/

&PARAM_ATMOS_VARS
 ATMOS_RESTART_OUTPUT         = .true.,
 ATMOS_RESTART_OUT_BASENAME   = "init_coldbubble",
/

&PARAM_MKINIT
 MKINIT_initname = "COLDBUBBLE",
/

&PARAM_MKINIT_COLDBUBBLE
 BBL_THETA = -5.D0,
 BBL_CZ = 6.D2,
 BBL_CX = 6.D2,
 BBL_CY = 6.D2,
 BBL_RZ = 2.D2,
 BBL_RX = 2.D2,
 BBL_RY = 2.D2,
/

End_of_SYSIN
########################################################################

# run
echo "job ${RUNNAME} started at " `date`
mpiexec.hydra -hostfile $TMPDIR/machines -n 1 -perhost 1 $BIN/$EXE ${EXE}.cnf
echo "job ${RUNNAME} end     at " `date`

exit
