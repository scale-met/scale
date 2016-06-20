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
INITNAME    = scale-rm_init
BINNAME     = scale-rm
PPCONF      = pp.conf
INITCONF    = init.conf
RUNCONF     = run.conf
TPROC       = 16


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
