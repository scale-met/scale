################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../../..)
TESTDIR     = ../../../..
DATADIR     = ${TOPDIR}/data
DATADIR_RM  = ${TESTDIR}/data

ifndef PREPDIR
   PREPDIR = ../pp
endif

ifndef INITDIR
   INITDIR = ../init
endif

# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
RUNCONF     = run.launch.conf
TPROC       = 6

# required data (parameters,distributed files)
DATPARAM    = " [${DATADIR}/rad/PARAG.29           PARAG.29         ] \
                [${DATADIR}/rad/PARAPC.29          PARAPC.29        ] \
                [${DATADIR}/rad/VARDATA.RM29       VARDATA.RM29     ] \
                [${DATADIR}/rad/cira.nc            cira.nc          ] \
                [${DATADIR}/rad/MIPAS/day.atm      day.atm          ] \
                [${DATADIR}/rad/MIPAS/equ.atm      equ.atm          ] \
                [${DATADIR}/rad/MIPAS/sum.atm      sum.atm          ] \
                [${DATADIR}/rad/MIPAS/win.atm      win.atm          ] \
                [${DATADIR}/land/param.bucket.conf param.bucket.conf] "
DATDISTS    = " [1 $(PREPDIR)/topo_d01    topo_d01   ] \
                [2 $(PREPDIR)/topo_d02    topo_d02   ] \
                [1 $(PREPDIR)/landuse_d01 landuse_d01] \
                [2 $(PREPDIR)/landuse_d02 landuse_d02] \
                [1 $(INITDIR)/boundary_d01 boundary_d01] \
                [2 $(INITDIR)/boundary_d02 boundary_d02] \
                [1 $(INITDIR)/init_d01_20070714-180000.000 init_d01_20070714-180000.000] \
                [2 $(INITDIR)/init_d02_20070714-180000.000 init_d02_20070714-180000.000] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
