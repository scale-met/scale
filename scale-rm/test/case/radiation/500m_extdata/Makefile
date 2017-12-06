################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..


# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
INITCONF    = init.conf
RUNCONF     = run.conf
TPROC       = 1

# required data (parameters,distributed files)
DATPARAM    = " [${TESTDIR}/data/rad/PARAG.29      PARAG.29    ] \
                [${TESTDIR}/data/rad/PARAPC.29     PARAPC.29   ] \
                [${TESTDIR}/data/rad/VARDATA.RM29  VARDATA.RM29] \
                [${TESTDIR}/data/rad/cira.nc       cira.nc     ] \
                [${TESTDIR}/data/rad/MIPAS/day.atm day.atm     ] \
                [${TESTDIR}/data/rad/MIPAS/equ.atm equ.atm     ] \
                [${TESTDIR}/data/rad/MIPAS/sum.atm sum.atm     ] \
                [${TESTDIR}/data/rad/MIPAS/win.atm win.atm     ] "
DATDISTS    = " [${TPROC} ../500m/rad_offline_19000321-000000.000 rad_offline_19000321-000000.000] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
