################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../../..)
TESTDIR     = ../../../..

ifndef PREPDIR
   PREPDIR = ../pp
endif

# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
INITCONF    = init.launch.conf
TPROC       = 3

# required data (parameters,distributed files)
DATPARAM    = " [${TESTDIR}/data/land/param.bucket.conf param.bucket.conf ] \
                [${TESTDIR}/tutorial/real/data/namelist.grads_boundary.FNL.2005053112-2016051106 namelist.grads_boundary.FNL.2005053112-2016051106] "
DATDISTS    = " [1 $(PREPDIR)/topo_d01    topo_d01   ] \
                [2 $(PREPDIR)/topo_d02    topo_d02   ] \
                [1 $(PREPDIR)/landuse_d01 landuse_d01] \
                [2 $(PREPDIR)/landuse_d02 landuse_d02] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

# special execute
run: jobshell
	bash link-input.sh
	$(JOBSUB) run.sh
