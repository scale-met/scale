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
INITCONF    = init_launch.conf

TPROC       = 8

# required data (parameters,distributed files)
DATPARAM    = " [${DATADIR}/land/param.bucket.conf param.bucket.conf ] \
                [${TESTDIR}/tutorial/real/data/namelist.grads_boundary.FNL.2005053112-2015011400 namelist.grads_boundary.FNL.2005053112-2015011400] "
DATDISTS    = " [4 $(PREPDIR)/topo_d01    topo_d01   ] \
                [4 $(PREPDIR)/topo_d02    topo_d02   ] \
                [4 $(PREPDIR)/landuse_d01 landuse_d01] \
                [4 $(PREPDIR)/landuse_d02 landuse_d02] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

# special execute
run: jobshell
	bash link-input.sh
	$(JOBSUB) run.sh
