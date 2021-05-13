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
TPROC       = 8

# required data (parameters,distributed files)
DATPARAM    = " [../../mkllmap/gl05rl01pe08_t42 llmap] "
DATDISTS    = " [${TPROC} ../../mkiimap/gl04rl00pe02_gl05rl01pe08/iimap_GL04RL00pe02_GL05RL01pe08 iimap_GL04RL00pe02_GL05RL01pe08] \
                [2        ../../mkllmap/gl04rl00pe02_i90j45/sample    sample_GL04RL00] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
