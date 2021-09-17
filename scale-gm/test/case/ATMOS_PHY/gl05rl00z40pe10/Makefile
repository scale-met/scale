################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..
TESTDIR_RM  = ${TOPDIR}/scale-rm/test/

# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
RUNCONF     = run.conf
ICO2LLCONF  = ico2ll.conf
TPROC       = 10

# required data (parameters,distributed files)
DATPARAM    = " [${TESTDIR}/data/grid/vgrid/vgrid40_600m_24km.dat vgrid40_600m_24km.dat] \
                [$(SCALE_DB)/grid/llmap/gl05rl00pe10_i180j90 llmap]                      \
		[${TESTDIR_RM}/data/rad/PARAG.29      PARAG.29    ]                      \
                [${TESTDIR_RM}/data/rad/PARAPC.29     PARAPC.29   ]                      \
                [${TESTDIR_RM}/data/rad/VARDATA.RM29  VARDATA.RM29]                      \
                [${TESTDIR_RM}/data/rad/cira.nc       cira.nc     ]                      \
                [${TESTDIR_RM}/data/rad/MIPAS/day.atm day.atm     ]                      \
                [${TESTDIR_RM}/data/rad/MIPAS/equ.atm equ.atm     ]                      \
                [${TESTDIR_RM}/data/rad/MIPAS/sum.atm sum.atm     ]                      \
                [${TESTDIR_RM}/data/rad/MIPAS/win.atm win.atm     ]                      \
"
DATDISTS    = " [${TPROC} $(SCALE_DB)/grid/boundary/gl05rl00pe10/boundary_GL05RL00 boundary_GL05RL00] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
