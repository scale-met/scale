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

TPROC       = 400

# required data (parameters,distributed files)
DATPARAM    = " [${DATADIR}/land/param.bucket.conf param.bucket.conf ] \
                [${TESTDIR}/tutorial/real/data/namelist.grads_boundary.FNL.1999080100-2005053106 namelist.grads_boundary.FNL.1999080100-2005053106] "
DATDISTS    = " [${TPROC} $(PREPDIR)/topo    topo   ] \
                [${TPROC} $(PREPDIR)/landuse landuse] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

# special execute
run: jobshell
	bash link-input.sh
	$(JOBSUB) run.sh
