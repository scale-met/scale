################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..


# user-defined source files
CODE_DIR    = ../code
ORG_SRCS    = mod_user.F90 mod_mkinit.F90  mod_atmos_phy_mp_driver.F90

# parameters for run
INITCONF    = init.conf
RUNCONF     = run.conf
TPROC       = 6

# required data (look-up table for charge split)
DATPARAM    = " [${TESTDIR}/data/lightning/LUT_TK1978_v.txt LUT_TK1978_v.txt] "
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
