################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..


# user-defined source files
CODE_DIR    = ..
ORG_SRCS    = mod_user.F90

# parameters for run
INITCONF    = init.conf
RUNCONF     = run.conf
TPROC       = 8

# required data (parameters,distributed files)
DATPARAM    = " [${TESTDIR}/data/input_sounding/twpice.txt   twpice.txt         ] \
                [${TESTDIR}/data/forcing/forcing_twpice.l120 forcing_twpice.l120] "
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
