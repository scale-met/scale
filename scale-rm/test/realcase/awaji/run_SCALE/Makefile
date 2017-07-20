################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

ifndef PREPDIR
   PREPDIR = ../pp
endif

ifndef INITDIR
   INITDIR = ../init_SCALE
endif

# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
RUNCONF     = run.conf
TPROC       = 16

# required data (parameters,distributed files)
DATPARAM    = " [${TESTDIR}/data/rad/PARAG.29           PARAG.29         ] \
                [${TESTDIR}/data/rad/PARAPC.29          PARAPC.29        ] \
                [${TESTDIR}/data/rad/VARDATA.RM29       VARDATA.RM29     ] \
                [${TESTDIR}/data/rad/cira.nc            cira.nc          ] \
                [${TESTDIR}/data/rad/MIPAS/day.atm      day.atm          ] \
                [${TESTDIR}/data/rad/MIPAS/equ.atm      equ.atm          ] \
                [${TESTDIR}/data/rad/MIPAS/sum.atm      sum.atm          ] \
                [${TESTDIR}/data/rad/MIPAS/win.atm      win.atm          ] \
                [${TESTDIR}/data/land/param.bucket.conf param.bucket.conf] "
DATDISTS    = " [${TPROC} $(PREPDIR)/topo     topo    ] \
                [${TPROC} $(PREPDIR)/landuse  landuse ] \
                [${TPROC} $(INITDIR)/boundary boundary] \
                [${TPROC} $(INITDIR)/init_19990301-060000.000 init_19990301-060000.000] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
