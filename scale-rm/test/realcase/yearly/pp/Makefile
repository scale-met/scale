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
PPCONF      = pp.conf

TPROC       = 400

# required data (parameters,distributed files)
DATPARAM    = " [$(SCALE_DB)/topo/GTOPO30/Products   input_topo   ] \
                [$(SCALE_DB)/landuse/GLCCv2/Products input_landuse] "
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
