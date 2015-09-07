################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..
SYSDEP_DIR  = $(abspath ..)

# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
PPNAME      = NONE
INITNAME    = scale-les_init
BINNAME     = scale-les
PPCONF      = NONE
INITCONF    = init.conf
RUNCONF     = run.conf
TPROC       = 24

# required data (parameters,distributed files)
DATDIR      =
DATPARAM    =
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
