#! /bin/bash -x
#
# for Local machine (MacOSX 8core+ifort+MPICH2)
#
export HMDIR=~/GCMresults/sol/latest
export BIN=~/Dropbox/Inbox/scale3/bin/${SCALE_SYS}
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
/usr/local/mpich213/bin/mpiexec -np 1 -f /Users/yashiro/libs/mpilib/machines_local $BIN/$EXE ${EXE}.cnf > STDOUT 2>&1
echo "job ${RUNNAME} end     at " `date`

exit
