################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

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
DATPARAM    = " [${TESTDIR}/data/land/param.bucket.conf param.bucket.conf ] \
                [${TESTDIR}/tutorial/real/data/namelist.grads_boundary.FNL.2005053112-2016051106 namelist.grads_boundary.FNL.2005053112-2016051106] "
DATDISTS    = " [$(PREPDIR)/topo    topo   ] \
                [$(PREPDIR)/landuse landuse] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

# special execute
run: jobshell
	bash link-input.sh
	$(JOBSUB) run.sh