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
PPNAME      = scale-rm_pp
INITNAME    = NONE
BINNAME     = NONE
PPCONF      = pp.conf
INITCONF    = NONE
RUNCONF     = NONE
TPROC       = 9

# required data (parameters,distributed files)
#DATDIR      = $(abspath /data2/scale_database/topo/DEM50M/Products)
DATDIR      = $(abspath /data2/scale_database/topo/GTOPO30/Products)
DATPARAM    =
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
