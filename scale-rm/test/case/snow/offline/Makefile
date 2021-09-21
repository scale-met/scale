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
ORG_SRCS    = mod_user.F90

# parameters for run
INITCONF    = init.conf
RUNCONF     = run.conf
TPROC       = 1

# required data (parameters,distributed files)
DATPARAM    = " [${DATADIR}/land/param.bucket.conf param.bucket.conf] "
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
