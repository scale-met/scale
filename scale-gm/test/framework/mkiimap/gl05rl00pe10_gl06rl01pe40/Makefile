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
PPNAME      = gm_mkiimap$(POSTFIX)
PPCONF      = pp.conf
ICO2LLNAME  = gm_fio_ico2ll
ICO2LLCONF  = ico2ll_iimap.conf
TPROC       = 40

# required data (parameters,distributed files)
DATPARAM    = " [../../mkllmap/gl06rl01pe40_i360j180 llmap] "
DATDISTS    = " [${TPROC} ../../mkhgrid/gl06rl01pe40/boundary_GL06RL01 boundary_GL06RL01] \
                [10       ../../mkhgrid/gl05rl00pe10/boundary_GL05RL00 boundary_GL05RL00] "


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
