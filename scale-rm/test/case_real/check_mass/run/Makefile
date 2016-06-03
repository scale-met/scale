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
INITNAME    = NONE
BINNAME     = scale-rm
PPCONF      = NONE
INITCONF    = NONE
RUNCONF     = run.conf
TPROC       = 4

# required data (parameters,distributed files)
DATDIR      = ../../../data
DATPARAM    = rad/PARAG.29,rad/PARAPC.29,rad/VARDATA.RM29,rad/cira.nc,rad/MIPAS/day.atm,rad/MIPAS/equ.atm,rad/MIPAS/sum.atm,rad/MIPAS/win.atm,land/param.bucket.conf
DATDISTS    = $(PREPDIR)/topo/topo.d01,$(PREPDIR)/landuse/landuse.d01,$(PREPDIR)/init/boundary.d01,$(PREPDIR)/init/init.d01_00011685600.000,$(PREPDIR)/topo/topo.d02,$(PREPDIR)/landuse/landuse.d02,$(PREPDIR)/init/init.d02_00011685600.000


export CHECK_MASS = T

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

run:
	$(JOBSUB) run.sh
