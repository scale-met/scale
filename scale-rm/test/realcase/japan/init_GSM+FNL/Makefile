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

ifndef PREPDIR
   PREPDIR = ../pp
endif

# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
INITCONF    = init.conf

TPROC       = 16

# required data (parameters,distributed files)
DATPARAM    = " [${DATADIR}/land/param.bucket.conf param.bucket.conf ] \
                [${PWD}/namelist.grads_boundary.GSM+FNL namelist.grads_boundary.GSM+FNL] "
DATDISTS    = " [${TPROC} $(PREPDIR)/topo    topo   ] \
                [${TPROC} $(PREPDIR)/landuse landuse] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

# special execute
run: jobshell
	bash link-input.sh
	$(JOBSUB) run.sh
