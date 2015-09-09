################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

# for build
index      = 200x8x8
precision  = real8
dynamics   = heve
tracer     = sdm
sfcflx     = const
turbulence = smg
radiation  = dymmy

# for run
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 6

# srcfiles
CODE_DIR = $(PWD)/code
ORG_SRCS = mod_user.f90 mod_mkinit.f90


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
