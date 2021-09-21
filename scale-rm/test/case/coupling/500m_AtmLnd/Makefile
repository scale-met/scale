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
DATPARAM    = " [${DATADIR}/rad/PARAG.29           PARAG.29         ] \
                [${DATADIR}/rad/PARAPC.29          PARAPC.29        ] \
                [${DATADIR}/rad/VARDATA.RM29       VARDATA.RM29     ] \
                [${DATADIR}/rad/cira.nc            cira.nc          ] \
                [${DATADIR}/rad/MIPAS/day.atm      day.atm          ] \
                [${DATADIR}/rad/MIPAS/equ.atm      equ.atm          ] \
                [${DATADIR}/rad/MIPAS/sum.atm      sum.atm          ] \
                [${DATADIR}/rad/MIPAS/win.atm      win.atm          ] \
                [${DATADIR}/land/param.bucket.conf param.bucket.conf] "
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
