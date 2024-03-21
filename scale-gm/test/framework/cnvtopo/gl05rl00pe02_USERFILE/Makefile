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
PPNAME      = gm_ll2ico$(POSTFIX)
PPCONF      = pp.conf
ICO2LLCONF  = ico2ll.conf
TPROC       = 2

# required data (parameters,distributed files)
DATPARAM    = " [$(SCALE_DB)/topo/userfile input_topo] \
                [$(SCALE_DB)/grid/llmap/gl05rl00pe02_i180j90 llmap] "
DATDISTS    = " [${TPROC} $(SCALE_DB)/grid/boundary/gl05rl00pe02/boundary_GL05RL00 boundary_GL05RL00] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
