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
PPNAME      = scale-les_pp
INITNAME    = NONE
BINNAME     = NONE
PPCONF      = pp.conf
INITCONF    = NONE
RUNCONF     = NONE
TPROC       = 36

# required data (parameters,distributed files)
DATDIR      = $(abspath /data2/scale_database/landuse/LU100M/Products)
DATPARAM    =
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include ../../Makefile.pp.common
