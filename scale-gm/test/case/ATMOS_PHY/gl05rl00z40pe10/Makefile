################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..
DATADIR     = ${TOPDIR}/data
DATADIR_GM  = ${TESTDIR}/data

# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
RUNCONF     = run.conf
ICO2LLCONF  = ico2ll.conf
TPROC       = 10

# required data (parameters,distributed files)
DATPARAM    = " [${DATADIR_GM}/grid/vgrid/vgrid40_600m_24km.dat vgrid40_600m_24km.dat] \
		[${DATADIR}/rad/PARAG.29      PARAG.29    ]                            \
                [${DATADIR}/rad/PARAPC.29     PARAPC.29   ]                            \
                [${DATADIR}/rad/VARDATA.RM29  VARDATA.RM29]                            \
                [${DATADIR}/rad/cira.nc       cira.nc     ]                            \
                [${DATADIR}/rad/MIPAS/day.atm day.atm     ]                            \
                [${DATADIR}/rad/MIPAS/equ.atm equ.atm     ]                            \
                [${DATADIR}/rad/MIPAS/sum.atm sum.atm     ]                            \
                [${DATADIR}/rad/MIPAS/win.atm win.atm     ]                            \
                [$(SCALE_DB)/grid/llmap/gl05rl00pe10_i180j90 llmap]                    \
"
DATDISTS    = " [${TPROC} $(SCALE_DB)/grid/boundary/gl05rl00pe10/boundary_GL05RL00 boundary_GL05RL00] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
