################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..
DATADIR     = ${TOPDIR}/data
DATADIR_RM  = ${TESTDIR}/data


# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
INITCONF    = init.conf
RUNCONF     = run.conf
TPROC       = 1

# required data (parameters,distributed files)
DATPARAM    = " [${DATADIR}/rad/PARAG.29      PARAG.29    ] \
                [${DATADIR}/rad/PARAPC.29     PARAPC.29   ] \
                [${DATADIR}/rad/VARDATA.RM29  VARDATA.RM29] \
                [${DATADIR}/rad/cira.nc       cira.nc     ] \
                [${DATADIR}/rad/MIPAS/day.atm day.atm     ] \
                [${DATADIR}/rad/MIPAS/equ.atm equ.atm     ] \
                [${DATADIR}/rad/MIPAS/sum.atm sum.atm     ] \
                [${DATADIR}/rad/MIPAS/win.atm win.atm     ] "
DATDISTS    = " [${TPROC} ../500m/rad_offline_19000321-000000.000 rad_offline_00000] \
                [${TPROC} ../500m/rad_offline_19000321-033000.000 rad_offline_00001] \
                [${TPROC} ../500m/rad_offline_19000321-063000.000 rad_offline_00002] \
                [${TPROC} ../500m/rad_offline_19000321-093000.000 rad_offline_00003] \
                [${TPROC} ../500m/rad_offline_19000321-123000.000 rad_offline_00004] \
                [${TPROC} ../500m/rad_offline_19000321-153000.000 rad_offline_00005] \
                [${TPROC} ../500m/rad_offline_19000321-183000.000 rad_offline_00006] \
                [${TPROC} ../500m/rad_offline_19000321-213000.000 rad_offline_00007] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
