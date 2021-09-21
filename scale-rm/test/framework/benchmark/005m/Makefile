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
SYSDEP_DIR  = $(abspath ..)

# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
INITCONF    = init.conf
RUNCONF     = run.conf
TPROC       = 6

# required data (parameters,distributed files)
DATPARAM    = " [${DATADIR}/rad/PARAG.29         PARAG.29        ] \
                [${DATADIR}/rad/PARAPC.29        PARAPC.29       ] \
                [${DATADIR}/rad/VARDATA.RM29     VARDATA.RM29    ] \
                [${DATADIR}/rad/rad_o3_profs.txt rad_o3_profs.txt] "
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
