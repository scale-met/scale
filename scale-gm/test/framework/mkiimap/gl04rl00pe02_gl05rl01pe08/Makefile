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
TPROC       = 8

# required data (parameters,distributed files)
DATPARAM    = " [../../mkllmap/gl05rl01pe08_t42 llmap] "
DATDISTS    = " [${TPROC} ../../mkhgrid/gl05rl01pe08/boundary_GL05RL01 boundary_GL05RL01] \
                [2        ../../mkhgrid/gl04rl00pe02/boundary_GL04RL00 boundary_GL04RL00] "


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
