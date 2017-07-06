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
   INITDIR = ../init_FNL
endif

# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
RUNCONF     = runA_launch.conf,runB1_launch.conf,runB2_launch.conf
TPROC       = 8,8,8

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
DATDISTS    = " [4 $(PREPDIR)/topo_d01    topo_d01   ] \
                [4 $(PREPDIR)/topo_d02    topo_d02   ] \
                [4 $(PREPDIR)/landuse_d01 landuse_d01] \
                [4 $(PREPDIR)/landuse_d02 landuse_d02] \
                [4 $(INITDIR)/boundary_d01 boundary_d01] \
                [4 $(INITDIR)/init_d01_20100501-000000.000 init_d01_20100501-000000.000] \
                [4 $(INITDIR)/init_d02_20100501-000000.000 init_d02_20100501-000000.000] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
