################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../..)
TESTDIR     = ../..

sample     ?= R20kmDX500m

# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
INITCONF    = ./sample/init_$(sample).conf
RUNCONF     = ./sample/run_$(sample).conf
N2GCONF     = ./sample/net2g_$(sample).conf
TPROC      ?= 2

# required data (parameters,distributed files)
DATDIR      =
DATPARAM    =
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
