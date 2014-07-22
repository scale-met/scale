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
PPNAME      = NONE
INITNAME    = NONE
BINNAME     = scale-les
PPCONF      = NONE
INITCONF    = NONE
RUNCONF     = run.conf
TPROC       = 24

# required data (parameters,distributed files)
DATDIR      = ../../../data
DATPARAM    = rad/PARAG.29,rad/PARAPC.29,rad/VARDATA.RM29,rad/cira.nc,rad/MIPAS/day.atm,rad/MIPAS/equ.atm,rad/MIPAS/sum.atm,rad/MIPAS/win.atm
DATDISTS    = /data2/converted_data/check_awajishima/topo/topo,/data2/converted_data/check_awajishima/landuse/landuse,/data2/converted_data/check_awajishima/initial_WRF/boundary,/data2/converted_data/check_awajishima/init_00000000000.000



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
