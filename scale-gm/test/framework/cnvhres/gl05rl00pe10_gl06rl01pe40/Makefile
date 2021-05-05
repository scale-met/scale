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
PPNAME     = gm_fio_ico2ico$(POSTFIX)
PPCONF     = ico2ico.conf
ICO2LLNAME  = gm_fio_ico2ll$(POSTFIX)
ICO2LLCONF  = ico2ll.conf
TPROC       = 40

# required data (parameters,distributed files)
DATPARAM    = " [../../mkllmap/gl06rl01pe40_i360j180 llmap] "
DATDISTS    = " [${TPROC} ../../mkiimap/gl05rl00pe10_gl06rl01pe40/iimap_GL05RL00pe10_GL06RL01pe40 iimap_GL05RL00pe10_GL06RL01pe40] \
                [10       ../../mkllmap/gl05rl00pe10_i180j90/sample    sample_GL05RL00] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
