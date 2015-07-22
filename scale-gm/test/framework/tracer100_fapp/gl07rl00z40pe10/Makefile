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
glevel      = 7
rlevel      = 0
nmpi        = 10
zlayer      = 40
vgrid       = vgrid40_600m_24km.dat



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
