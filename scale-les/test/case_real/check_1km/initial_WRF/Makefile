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
INITNAME    = scale-les_init
BINNAME     = NONE
PPCONF      = NONE
INITCONF    = init.conf
RUNCONF     = NONE
TPROC       = 16

# required data (parameters,distributed files)
DATDIR      = $(abspath ../../../../../../../input_database/WRF_output/kenshi_gtopo30_case1)
DATPARAM    =
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include ../../Makefile.pp.common
