################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..


# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
INITCONF    = init.conf
RUNCONF     = run.conf
TPROC       = 1

# required data (parameters,distributed files)
DATDIR      = ../../../data
DATPARAM    = rad/PARAG.29,rad/PARAPC.29,rad/VARDATA.RM29,rad/cira.nc,rad/MIPAS/day.atm,rad/MIPAS/equ.atm,rad/MIPAS/sum.atm,rad/MIPAS/win.atm
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
