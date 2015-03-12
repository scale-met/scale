################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

ifndef PREPDIR
PREPDIR     = ..
endif


# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
PPNAME      = NONE
INITNAME    = scale-les_init
BINNAME     = NONE
PPCONF      = NONE
INITCONF    = init.conf
RUNCONF     = NONE
TPROC       = 4

# required data (parameters,distributed files)
DATDIR      = ../../../data
DATPARAM    = rad/PARAG.29,rad/PARAPC.29,rad/VARDATA.RM29,rad/cira.nc,rad/MIPAS/day.atm,rad/MIPAS/equ.atm,rad/MIPAS/sum.atm,rad/MIPAS/win.atm,land/param.bucket.conf
DATDISTS    = $(PREPDIR)/topo/topo.d01,$(PREPDIR)/topo/topo.d02,$(PREPDIR)/landuse/landuse.d01,$(PREPDIR)/landuse/landuse.d02


export CHECK_MASS = T

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

# special execute
run: jobshell
	$(JOBSUB) nicaminput-link.sh
	$(JOBSUB) run.sh
